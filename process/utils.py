#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 15:51:31 2019

@author: benjamin.garzon
"""
import nibabel as nb
from nibabel.gifti import gifti, giftiio
from mvpa2.suite import *
from mvpa2.measures.rsa import PDist
from scipy import stats
from sklearn.decomposition import PCA

def distinctiveness(x):
    return np.mean(x, axis=0)

def myarctanh(x):
    x[np.isnan(x)] = 0
    x[np.abs(x) >= 1] = 0
    
#    print("x", x)
#    print("atan", np.arctanh(x))
    return(np.arctanh(x))

def get_spread_chunk(x, metric, targets):
    RDM = squareform(x)
    if metric == 'correlation':
        RDM = myarctanh(1 - RDM) # transform to normalized correlation

    elements = np.unique(targets)
    N = len(elements)
    within_values = {}
    between_values = {} 
    
    for i, elementi in enumerate(elements):
        seli = targets == elementi 
        M = np.sum(seli)

        if M > 1:
            within_values[elementi] = []
            between_values[elementi] = [] 
        
            for j, elementj in enumerate(elements):
                selj = targets == elementj 
                if i != j:
                    if np.sum(selj) > 0:
                        X = RDM[seli, :][:, selj].ravel()
                        between_values[elementi].append(X)
                else:
                    X = RDM[seli, :][:, seli][np.tril_indices(M, k = -1)]
                    within_values[elementi].append(X)
    
    return(within_values, between_values)

# faster, only looks at within values
def get_within_spread_chunk(x, metric, targets):
    RDM = squareform(x)
    if metric == 'correlation':
        RDM = myarctanh(1 - RDM)

    elements = np.sort(np.unique(targets))
    N = len(elements)
    within_values = {}
    for i, elementi in enumerate(elements):
        seli = targets == elementi 
        M = np.sum(seli)

        if M > 1:
            within_values[elementi] = []        
            X = RDM[seli, :][:, seli][np.tril_indices(M, k = -1)]
            within_values[elementi].append(X)
            
    return(within_values)
    
def get_RDM_metric(x, metric, targets):
    RDM = squareform(x)
    if metric == 'correlation':
        RDM = myarctanh(1 - RDM)
    elements = np.unique(targets)
    
    N = len(elements)
    meanRDM = np.zeros((N, N)) 
    
    within_values = []
    between_values = []

    for i, elementi in enumerate(elements):
        seli = targets == elementi 
        for j, elementj in enumerate(elements):
            selj = targets == elementj 
            if i != j:
                meanRDM[i, j] = np.mean(RDM[seli, :][:, selj])
                if i > j:
                    between_values.append(RDM[seli, :][:, selj].ravel())
            else:
                M = np.sum(seli)
                meanRDM[i, i] = np.mean(RDM[seli, :][:, seli][np.tril_indices(M, k = -1)])
                within_values.append(RDM[seli, :][:, seli][np.tril_indices(M, k = -1)])
    within = np.nanmean(np.concatenate(within_values))
    between = np.nanmean(np.concatenate(between_values))
        
#    if metric == 'correlation':
#        meanRDM = 1 - np.tanh(meanRDM)
#        within = 1 - np.tanh(within)
#        between = 1 - np.tanh(between)
    
    return(meanRDM, within, between)    

def get_spread(x, *args):
    meanRDM, within, between = get_RDM_metric(x, args[0], args[1])
    if args[0] != 'euclidean':
        spread = between - within
    else:
        spread = between/within
    return(spread)

def get_within(x, *args):
    meanRDM, within, between = get_RDM_metric(x, args[0], args[1])
    return(within)

def get_between(x, *args):
    meanRDM, within, between = get_RDM_metric(x, args[0], args[1])
    return(between)
      

distinct_mapper = FxMapper('samples', distinctiveness, attrfx='merge')

def plot_mtx(mtx, labels, title, clim=(0, 2)):
    pl.figure()
    pl.imshow(mtx, interpolation='nearest')
    pl.xticks(range(len(mtx)), labels, rotation=-45)
    pl.yticks(range(len(mtx)), labels)
    pl.title(title)
    pl.clim(clim)
    pl.colorbar()


# deprecated :
    
def aggregate_matrix_normalize(RDM, chunks, targets):
    
    RSM = np.arctanh(1 - RDM) # as correlation and Fisher transform it
    
    elements = np.unique(targets)
    
    unique_chunks = np.unique(chunks) 
    N = len(elements)
    Nchunks = len(unique_chunks)
    meanRSM = np.zeros((N, N, Nchunks)) 

    within = np.zeros(Nchunks)
    between = np.zeros(Nchunks)
    
    for k, chunk in enumerate(unique_chunks):    
        within_values = []
        between_values = []

        for i, elementi in enumerate(elements):
            seli = np.logical_and(targets == elementi, chunks == chunk) 
            for j, elementj in enumerate(elements):
                selj = np.logical_and(targets == elementj, chunks == chunk) 
                if i != j:
                    meanRSM[i, j, k] = np.mean(RSM[seli, :][:, selj])
                    if i > j:
                        between_values.append(RSM[seli, :][:, selj].ravel())
                else:
                    M = np.sum(seli)
                    meanRSM[i, i, k] = np.mean(RSM[seli, :][:, seli][np.tril_indices(M, k = -1)])
                    within_values.append(RSM[seli, :][:, seli][np.tril_indices(M, k = -1)])
#        pl.figure(figsize = (9, 6))    
#        pl.plot(np.concatenate(within_values), 'r')
#        pl.plot(np.concatenate(between_values), 'b')
#        pl.title(chunk)
        within[k] = 1 - np.tanh(np.nanmean(np.concatenate(within_values)))
        between[k] = 1 - np.tanh(np.nanmean(np.concatenate(between_values)))
        
    score = np.mean(within/between) # across chunks    
    meanRDM = 1 - np.tanh(meanRSM)
    return(meanRDM, elements, score, within, between)


def aggregate_matrix(RDM, chunks, targets):
    
    elements = np.unique(targets)
    
    unique_chunks = np.unique(chunks) 
    N = len(elements)
    Nchunks = len(unique_chunks)
    meanRDM = np.zeros((N, N, Nchunks)) 

    within = np.zeros(Nchunks)
    between = np.zeros(Nchunks)
    
    for k, chunk in enumerate(unique_chunks):    
        within_values = []
        between_values = []

        for i, elementi in enumerate(elements):
            seli = np.logical_and(targets == elementi, chunks == chunk) 
            for j, elementj in enumerate(elements):
                selj = np.logical_and(targets == elementj, chunks == chunk) 
                if i != j:
                    meanRDM[i, j, k] = np.mean(RDM[seli, :][:, selj])
                    if i > j:
                        between_values.append(RDM[seli, :][:, selj].ravel())
                else:
                    M = np.sum(seli)
                    meanRDM[i, i, k] = np.mean(RDM[seli, :][:, seli][np.tril_indices(M, k = -1)])
                    within_values.append(RDM[seli, :][:, seli][np.tril_indices(M, k = -1)])
                    
#        pl.figure(figsize = (9, 6))    
#        pl.plot(np.concatenate(within_values), 'r')
#        pl.plot(np.concatenate(between_values), 'b')
#        pl.title(chunk)
        within[k] = np.nanmean(np.concatenate(within_values))
        between[k] = np.nanmean(np.concatenate(between_values))
        
    score = np.mean(within/between) # across chunks    
    return(meanRDM, elements, score, within, between)


class flatten_mapper(Mapper):

    is_trained = True
    
    def __init__(self, **kwargs):
        Mapper.__init__(self, **kwargs)
    
    def _forward_dataset(self, ds):
        mapped_ds = ds.copy()
        # standardize so that the two dimensions are comparable
#        ds0 = stats.zscore(ds.samples[:, :, 0], axis = None)
#        ds1 = stats.zscore(ds.samples[:, :, 1], axis = None)
        ds0 = ds.samples[:, :, 0]
        ds1 = ds.samples[:, :, 1]
        mapped_ds.samples = np.hstack((ds0, ds1))
        return(mapped_ds)


def get_spread_chunk_map(x, *args):
    return(get_spread_chunk(x, args[0], args[1], args[2]))


class PClassifier(PDist):

    def __init__(self, 
                 classifier,
                 mapper = None, 
                 NCOMPS = 10, 
                 filter_accuracy = False,
                 accuracy = 0, 
                 **kwargs):
        self.classifier = classifier
        self.mapper = mapper
        self.filter_accuracy = filter_accuracy     
        self.correct = accuracy == 1
        self.pca = PCA(n_components = NCOMPS, whiten = False)
        super(PDist, self).__init__(**kwargs)
        
    def _call(self, ds):

        mapped_ds = self.mapper(ds)
        data = mapped_ds.samples
        # remove constant columns
        constantcols = np.all(data[1:] == data[:-1], axis=0)
        data = data[:, ~constantcols]
        if np.sum(constantcols)>0:
            print(np.sum(constantcols))
        try:
            data_pca = self.pca.fit_transform(data)
            
            if self.filter_accuracy:
                data_pca = data_pca[self.correct]
                chunks = ds.chunks[self.correct]
                targets = ds.targets[self.correct]
            else:
                chunks = ds.chunks
                targets = ds.targets
   
            data_red = Dataset(data_pca, 
                               sa={'chunks': chunks, 'targets': targets})
            
            out = self.classifier(data_red)

        except LinAlgError:
            print("PCA did not converge!")
            out = Dataset(np.zeros((len(np.unique(ds.targets[self.correct])), 1)))           
            
#        print(out.samples)          
        return out



class Pspread(PDist):

    def __init__(self, 
                 mapper = None, 
                 NCOMPS = 10, 
                 filter_accuracy = False,
                 accuracy = 0,
                 **kwargs):

        self.mapper = mapper
        self.filter_accuracy = filter_accuracy     
        self.correct = accuracy == 1
        self.pca = PCA(n_components = NCOMPS, whiten = False)
        super(Pspread, self).__init__(**kwargs)

    def _call(self, ds):

        mapped_ds = self.mapper(ds)
        data = mapped_ds.samples
        # remove constant columns
        constantcols = np.all(data[1:] == data[:-1], axis=0)
        data = data[:, ~constantcols]
        if np.sum(constantcols)>0:
            print(np.sum(constantcols))
        try:
            data_pca = self.pca.fit_transform(data)
            
            if self.filter_accuracy:
                data_pca = data_pca[self.correct]
                chunks = ds.chunks[self.correct]
                targets = ds.targets[self.correct]
            else:
                chunks = ds.chunks
                targets = ds.targets
#                print(data_pca.shape)
    
            unique_chunks = np.unique(chunks) 
            Nchunks = len(unique_chunks)
            within_mean = {}
            between_mean = {}
    
            for k, chunk in enumerate(unique_chunks):
                ds_chunk = data_pca[chunks == chunk]
                dsm = pdist(ds_chunk, 
                        metric=self.params.pairwise_metric)
                within, between = get_spread_chunk(dsm, 
                                   self.params.pairwise_metric, 
                                   targets[chunks == chunk])
                
                if k == 0:
                    within_values = {}
                    between_values = {}
                    
                for key in within.keys():
                    if key in within_values.keys():
                        within_values[key] = within_values[key] + within[key]
                    else:
                        within_values[key] = within[key]

                            
                for key in between.keys():
                    if key in between_values.keys():
                        between_values[key] = between_values[key] + between[key]
                    else:
                        between_values[key] = between[key]
                    
            for key in within_values.keys():
                within_mean[key] = np.nanmean(np.concatenate(within_values[key]))                        
                between_mean[key] = np.nanmean(np.concatenate(between_values[key]))
                                     
            if self.params.pairwise_metric == 'correlation':
                spread = np.mean(
                        (1 - np.tanh(np.array(between_mean.values())))/
                        (1 - np.tanh(np.array(within_mean.values())))
                        )
                        
            else:
                spread = np.mean(np.array(between_mean.values())/
                                         np.array(within_mean.values()))
        except LinAlgError:
            print("PCA did not converge!")
 #           print(data)
            spread = 0
        except ValueError:
            print("Value error")
 #           print(data)
            spread = 0
            
        out = Dataset(np.array((spread,)))
        return out


class Pwithin_spread(PDist):

    def __init__(self, 
                 seq_train,
                 mapper = None, 
                 NCOMPS = 10, 
                 filter_accuracy = False,
                 accuracy = 0,                 
                 **kwargs):

        self.mapper = mapper
        self.filter_accuracy = filter_accuracy     
        self.correct = accuracy == 1
        self.seq_train = seq_train
        self.pca = PCA(n_components = NCOMPS, whiten = False)
        super(Pwithin_spread, self).__init__(**kwargs)

    def _call(self, ds):

        mapped_ds = self.mapper(ds)
        data = mapped_ds.samples
        
        # remove constant columns
        constantcols = np.all(data[1:] == data[:-1], axis=0)
        data = data[:, ~constantcols]
        if np.sum(constantcols)>0:
            print(np.sum(constantcols))
        try:
            data_pca = self.pca.fit_transform(data)
            
            if self.filter_accuracy:
                data_pca = data_pca[self.correct]
                chunks = ds.chunks[self.correct]
                targets = ds.targets[self.correct]
            else:
                chunks = ds.chunks
                targets = ds.targets
#                print(data_pca.shape)
    
            unique_chunks = np.unique(chunks) 
            Nchunks = len(unique_chunks)
            within_mean = {}
    
            for k, chunk in enumerate(unique_chunks):
                ds_chunk = data_pca[chunks == chunk]
                dsm = pdist(ds_chunk, 
                        metric=self.params.pairwise_metric)
                within = get_within_spread_chunk(dsm, 
                                   self.params.pairwise_metric, 
                                   targets[chunks == chunk])
                if k == 0:
                    within_values = {}
                    
                for key in within.keys():
                    if key in within_values.keys():
                        within_values[key] = within_values[key] + within[key]
                    else:
                        within_values[key] = within[key]
                        
            within_trained = []
            within_untrained = []            
            
            for key in within_values.keys():
                within_mean[key] = np.nanmean(np.concatenate(within_values[key]))                        
                if self.seq_train[key] == 'untrained':
                    within_untrained.append(within_mean[key])
                else:
                    within_trained.append(within_mean[key])
                                     
            if self.params.pairwise_metric == 'correlation':
                within_spread = (1 - np.tanh(np.nanmean(within_untrained)))/ \
                (1 - np.tanh(np.nanmean(within_trained)))
                        
            else: # not ready
                within_spread = 0

        except LinAlgError:
            print("PCA did not converge!")
 #           print(data)
            spread = 0
        except ValueError:
            print("Value error")
 #           print(data)
            spread = 0
            
        out = Dataset(np.array((within_spread,)))
        return out
#        except LinAlgError:
#            print("PCA did not converge!")
# #           print(data)
#            within_spread = np.zeros((len(within_values.keys()), ))
#        except ValueError:
#            print("Value error")
# #           print(data)
#            within_spread = np.zeros((len(within_values.keys()), ))
#        out = Dataset(np.array((within_spread, )).T)
#        return out



class PDistMulti(PDist):

    def __init__(self, mapper = None, **kwargs):
        self.mapper = mapper
        super(PDistMulti, self).__init__(**kwargs)

    def _call(self, ds):

        mapped_ds = self.mapper(ds)
#        print(ds.shape, mapped_ds.shape)
#        print(mapped_ds)
        data = mapped_ds.samples
        # get dsm
        dsm = pdist(data, metric=self.params.pairwise_metric)

        # add some attributes
        out = Dataset(dsm, 
                      sa=dict(pairs=list(combinations(range(len(mapped_ds)), 2))))
        return out
#        out = super(PDistMulti, self)._call(self, mapped_ds)
#        print(out)
#        return(out)

class SearchlightPreproc(Searchlight):
    """Implements a searchlight with an additional preprocessing step so that 
    multidimensional features can be run.
    """
    
    def __init__(self, datameasure, queryengine, add_center_fa=False,
             results_postproc_fx=None,
             results_backend='native',
             results_fx=None,
             tmp_prefix='tmpsl',
             nblocks=None,
             preallocate_output=False,
             mapper = None, 
             **kwargs):
        self.mapper = mapper
        
        Searchlight.__init__(self, datameasure, 
                             queryengine, 
                             add_center_fa,
                             results_postproc_fx,
                             results_backend,
                             results_fx,
                             tmp_prefix,
                             nblocks,
                             preallocate_output,
                             **kwargs)
        
    def _call(self, dataset):
        
        mapped_dataset = self.mapper()(dataset)
        Searchlight._call(self, mapped_dataset)

def map2gifti2(ds, filename=None, encoding='GIFTI_ENCODING_B64GZ',
              surface=None, vertices=None):
    """Maps data(sets) into a GiftiImage, and optionally saves it to disc.
    Parameters
    ----------
    ds : AttrDataset or numpy.ndarray
      The data to be mapepd
    filename : basestring or None, optional
      Filename to which the GiftiImage is stored
    encoding : "ASCII" or "Base64Binary" or "GZipBase64Binary", optional
      Encoding format of data
    surface : mvpa2.surf.nibabel.surf.Surface or str, optional
      Optional anatomical Surface object, or filename of anatomical surface
      file, to be stored together with the data. This should allow
      FreeSurfer's mris_convert to read files written by this function
    Returns
    -------
    img : GiftiImage
      dataset contents represented in GiftiImage
    """

    darrays = []

    if isinstance(ds, np.ndarray):
        samples = ds
    elif isinstance(ds, AttrDataset):
        samples = ds.samples
        _warn_if_fmri_dataset(ds)
    else:
        raise TypeError('first argument must be AttrDataset or numpy.ndarray')

    [nsamples, nfeatures] = samples.shape

    def _get_attribute_value(ds, attr_name, keys_):
        if isinstance(ds, np.ndarray):
            # no attributes
            return None

        attr_collection = ds.__dict__.get(attr_name)

        if isinstance(keys_, basestring):
            keys_ = (keys_,)

        for key in keys_:
            if key in attr_collection:
                return attr_collection[key].value
        return None

    def _build_array(data, intent, encoding=encoding):
        is_integer = intent == 'NIFTI_INTENT_NODE_INDEX'
        dtype = np.int32 if is_integer else np.float32
        arr = gifti.GiftiDataArray(data.astype(dtype), intent,
                                              encoding=encoding)
        # Setting the coordsys argument the constructor would set the matrix
        # to the 4x4 identity matrix, which is not desired. Instead the
        # coordsys is explicitly set to None afterwards
        arr.coordsys = None

        return arr

    node_indices_labels = ('node_indices', 'center_ids', 'ids', 'roi_ids')
    node_indices = _get_attribute_value(ds, 'fa', node_indices_labels)
    
    if vertices is not None:
        values = np.zeros((samples.shape[0], vertices))
        for i, sample in enumerate(samples):
            values[i, node_indices] = sample
        darray = _build_array(np.arange(vertices), 'NIFTI_INTENT_NODE_INDEX')
        darrays.append(darray)
        samples = values

    else:
        if node_indices is not None:
            darray = _build_array(node_indices, 'NIFTI_INTENT_NODE_INDEX')
            darrays.append(darray)

    intents = _get_attribute_value(ds, 'sa', 'intents')
    for i, sample in enumerate(samples):
        intent = 'NIFTI_INTENT_NONE' if intents is None else intents[i]
        darray = _build_array(sample, intent)
        darrays.append(darray)

    # if there is a surface, add it
    if surface is not None:
        surface_object = surf_from_any(surface, )
        anat_image = anat_surf_to_gifti_image(surface_object, add_indices=False)

        for darray in anat_image.darrays:
            darrays.append(darray)

    image = gifti.GiftiImage(darrays=darrays)

    if filename is not None:
        nb.save(image, filename)


def _warn_if_fmri_dataset(ds):
    assert (isinstance(ds, AttrDataset))

    fmri_fields = set(('imgaffine', 'imgtype', 'imghdr'))

    ds_fmri_fields = set.intersection(set(ds.a.keys()), fmri_fields)

    if len(ds_fmri_fields) > 0:
        warning('dataset attribute .a has fields %s, which suggest it is an '
                'volumetric dataset. Converting this dataset to GIFTI '
                'format will most likely result in unvisualiable '
                '(and potentially, un-analysable) data. Consider using '
'map2nifti instead' % (', '.join(ds_fmri_fields)))
        
        
def mahalanobis_distance_regularized(x, y=None, w=None):
    """Calculate Mahalanobis distance of the pairs of points.
    Parameters
    ----------
    `x`
      first list of points. Rows are samples, columns are
      features.
    `y`
      second list of points (optional)
    `w` : np.ndarray
      optional inverse covariance matrix between the points. It is
      computed if not given
    Inverse covariance matrix can be calculated with the following
      w = np.linalg.solve(np.cov(x.T), np.identity(x.shape[1]))
    or
      w = np.linalg.inv(np.cov(x.T))
    """
    # see if pairwise between two matrices or just within a single matrix
    if y is None:
        # pairwise distances of single matrix
        # calculate the inverse correlation matrix if necessary
        if w is None:
            w = np.linalg.inv(np.cov(x.T))

        # get some shapes of the data
        mx, nx = x.shape
        #mw, nw = w.shape

        # allocate for the matrix to fill
        d = np.zeros((mx, mx), dtype=np.float32)
        for i in range(mx-1):
            # get the current row to compare
            xi = x[i, :]
            # replicate the row
            xi = xi[np.newaxis, :].repeat(mx-i-1, axis=0)
            # take the distance between all the matrices
            dc = x[i+1:mx, :] - xi
            # scale the distance by the correlation
            d[i+1:mx, i] = np.real(np.sum((np.inner(dc, w) * np.conj(dc)), 1))
            # fill the other direction of the matrix
            d[i, i+1:mx] = d[i+1:mx, i].T
    else:
        # is between two matrixes
        # calculate the inverse correlation matrix if necessary
        if w is None:
            # calculate over all points
            w = np.linalg.inv(np.cov(np.concatenate((x, y)).T))

        # get some shapes of the data
        mx, nx = x.shape
        my, ny = y.shape

        # allocate for the matrix to fill
        d = np.zeros((mx, my), dtype=np.float32)

        # loop over shorter of two dimensions
        if mx <= my:
            # loop over the x patterns
            for i in range(mx):
                # get the current row to compare
                xi = x[i, :]
                # replicate the row
                xi = xi[np.newaxis, :].repeat(my, axis=0)
                # take the distance between all the matrices
                dc = xi - y
                # scale the distance by the correlation
                d[i, :] = np.real(np.sum((np.inner(dc, w) * np.conj(dc)), 1))
        else:
            # loop over the y patterns
            for j in range(my):
                # get the current row to compare
                yj = y[j, :]
                # replicate the row
                yj = yj[np.newaxis, :].repeat(mx, axis=0)
                # take the distance between all the matrices
                dc = x - yj
                # scale the distance by the correlation
                d[:, j] = np.real(np.sum((np.inner(dc, w) * np.conj(dc)), 1))

    # return the dist
    return np.sqrt(d)


class Crossnobis(Measure):
    """Compute cross-validated dissimiliarity matrix for samples in a dataset
    This `Measure` can be trained on part of the dataset (for example,
    a partition) and called on another partition. It can be used in
    cross-validation to generate cross-validated RSA.
    Returns flattened dissimilarity values.
    """
    pairwise_metric = Parameter('correlation', constraints='str',
            doc="""Distance metric to use for calculating pairwise vector distances for
            dissimilarity matrix (DSM).  See scipy.spatial.distance.cdist for
            all possible metrics.""")

    pairwise_metric_kwargs = Parameter({},
            doc="""kwargs dictionary passed to cdist. For example,
            if `pairwise_metric='mahalanobis'`, `pairwise_metric_kwargs`
            might contain the inverse of the covariance matrix.""")

    sattr = Parameter(['targets'],
            doc="""List of sample attributes whose unique values will be used to
            identify the samples groups. Typically your category labels or targets.""")

    def __init__(self, **kwargs):
        Measure.__init__(self, **kwargs)
        self._train_ds = None

    def _prepare_ds(self, ds):
        if self.params.sattr is not None:
            mgs = mean_group_sample(attrs=self.params.sattr)
            ds_ = mgs(ds)
        else:
            ds_ = ds.copy()
        return ds_

    def _train(self, ds):
        self._train_ds = self._prepare_ds(ds)
        self.is_trained = True

    def _call(self, ds):
        test_ds = self._prepare_ds(ds)
        if test_ds.nsamples != self._train_ds.nsamples:
            raise ValueError('Datasets should have same sample size for dissimilarity, '\
                             'nsamples for train: %d, test: %d'%(self._train_ds.nsamples,
                                                                 test_ds.nsamples))
        # Call actual distance metric
        distds = cdist(self._train_ds.samples, test_ds.samples,
                       metric=self.params.pairwise_metric,
                       **self.params.pairwise_metric_kwargs)
        # Make target pairs
        sa_dict = dict()
        for k in self._train_ds.sa:
            if k in test_ds.sa:
                sa_dict[k] = list(product(self._train_ds.sa.get(k).value,
                                                   test_ds.sa.get(k).value))

        distds = Dataset(samples=distds.ravel()[:, None], sa=sa_dict)
        return distds

