from collections import Counter
import cPickle
import imp
import logging
import re
import pandas as pd
from suds.client import Client

class DAVID(object) :
  """Interface to the DAVID gene set enrichment toolset API
  
  Args:
    username (str): username (usually email address) used to authenticate with
      the DAVID API
    ids (list of str): gene identifiers consistent with *idType* argument
    idType (str): The ID type of the identifiers in *id*, can be any of the
      accepted ID types listed on the DAVID site
    summary_n (int): the number of most common terms to report in the summary
      functions
    stringency (str): the stringency argument, one of [Lowest, Low, Medium, High,
      Highest]
    load_fn (str): path to file that was saved using the DAVID.save method to
      load a pre-run analysis

  Methods:
    getGeneClusterReport(stringency="medium")
    getGeneClusterSummary(n=4)
    getTermClusterReport(stringency="medium")
    getTermClusterSummary(n=4)

  """
  def __init__(self,username,
         ids=None,
         idType='ENSEMBL_GENE_ID',
         summary_n=4,
         stringency='medium',
         load_fn=None) :
    
    if load_fn is not None : # load the DAVID analysis from a previously saved file
      (self.gene_cluster,
       self.gene_summary,
       self.term_cluster,
       self.term_summary) = cPickle.load(open(load_fn))
    else :
      if ids is None :
        raise Exception("Must supply a list of IDs if load_fn is not provided")
        
      if len(ids) > 3000 :
        raise Exception("DAVID API does not accept queries with > 3000 identifiers")
      
      ids = ', '.join(ids)
        
      url = 'http://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl'
    
      self.client = client = Client(url)
      client.set_options(timeout=600)
      logging.getLogger('suds.client').setLevel(logging.WARNING)
      
      client.service.authenticate(username)
      
      # add list of IDs
      listName = 'make_up'
      listType = 0
      client.service.addList(ids,idType,listName,listType)
      
      #self.chart = self.getChartReport()
      self.gene_cluster = self.getGeneClusterReport(stringency=stringency)
      self.gene_summary = self.getGeneClusterSummary(summary_n)
      self.term_cluster = self.getTermClusterReport(stringency=stringency)
      self.term_summary = self.getTermClusterSummary(summary_n)
    
  def _map_stringency(self,stringency,args=None) :
    strng = stringency.lower()
    if args is not None :
      pass
    elif strng == 'highest' :
      args = 5, 5, 5, 0.5, 100
    elif strng == 'high' :
      args = 4, 5, 5, 0.5, 85
    elif strng == 'medium' :
      args = 4, 4, 4, 0.5, 50
    elif strng == 'low' :
      args = 4, 3, 3, 0.5, 35
    elif strng == 'lowest' :
      args = 3, 3, 3, 0.5, 20
    else :
      args = 3, 3, 3, 0.5, 50
    return args

  def getChartReport(self,counts=2,threshold=0.1) :
    return self.client.service.getChartReport(threshold,counts)
    
  def getGeneClusterSummary(self,n=4) :
    return self._summary_report(self.gene_cluster,n)

  def getGeneClusterReport(self,stringency="medium",args=None) :
    # args = overlap, initialSeed, finalSeed, linkage, kappa
    args = self._map_stringency(stringency,args)
    self._gene_result = result = self.client.service.getGeneClusterReport(*args)
    clusters = []
    print 'found %d results, processing'%len(result)
    for c in result:
      df = pd.DataFrame(dict((r.values[0],dict(r.geneObject)) for r in c.listRecords)).T
      df.cluster_name = c.name
      df.cluster_score = c.score
      clusters.append(df)
    return clusters

  def getTermClusterSummary(self,n=4) :
    return self._summary_report(self.term_cluster,n)

  def _summary_report(self,res,n=4) :
    sub_patts = ('GO:\d+~','IPR\d+:','SM\d+:','domain:',
           'repeat:',',',':','[()]')
    stop_patts = ('of','in','a','to','on','\d+','and')
    freq_table = []
    for c in res :
      cntr = Counter()
      for w in c.name.tolist() :
        for sub_patt in sub_patts :
          w = re.sub(sub_patt,'',w)
        scrub_w = w.lower()
        for w_i in scrub_w.split() :
          if all([re.match('^%s$'%re.escape(p),re.escape(w_i)) == None for p in stop_patts]) :
            orig_w_i = re.search(re.escape(w_i),re.escape(w),re.IGNORECASE)
            g = orig_w_i
            if g is None :
              cntr[w_i] += 1
            else :
              cntr[orig_w_i.group()] += 1
      if len(cntr) != 0 :
        freq_table.append((c.cluster_name,', '.join(zip(*cntr.most_common(n))[0]),len(c.index),c.cluster_score))
        #freq_table.append((c.cluster_name,', '.join(cntr.most_common(n)),len(c.index),c.cluster_score))

    return pd.DataFrame(freq_table,columns=['Name','Top 4 terms','# terms','score'])

  def getTermClusterReport(self,stringency="medium",args=None) :
    # args = overlap, initialSeed, finalSeed, linkage, kappa
    args = self._map_stringency(stringency,args)
    self._term_result = result = self.client.service.getTermClusterReport(*args)
    clusters = []
    for c in result:
      
      df = pd.DataFrame(dict((r.termName,dict(r)) for r in c.simpleChartRecords)).T
      df['name'] = df.index
      df.cluster_name = c.name
      df.cluster_score = c.score
      clusters.append(df)
    return clusters
  
  def save(self,fn) :
    cPickle.dump([self.gene_cluster,
            self.gene_summary,
            self.term_cluster,
            self.term_summary],open(fn,'w'))
