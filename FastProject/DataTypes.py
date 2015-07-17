# -*- coding: utf-8 -*-
"""
Created on Thu Mar 05 20:13:40 2015

@author: daved_000
"""

import numpy as np;
from scipy.spatial import distance;
import scipy.stats;

class ExpressionData(np.ndarray):
    
    def __new__(subtype, data, row_labels=[], col_labels=[]):       
        
        obj = np.asarray(data).view(subtype);
        
        if(len(row_labels) == 0):
            obj.row_labels = [str(i) for i in range(data.shape[0])];
        elif(len(row_labels) != data.shape[0]):
            raise ValueError("Number of row labels does not equal number of rows in data");
        else:
            obj.row_labels = row_labels;
        
        if(len(col_labels) == 0):
            obj.col_labels = [str(i) for i in range(data.shape[1])];
        elif(len(col_labels) != data.shape[1]):
            raise ValueError("Number of column labels does not equal number of columns in data");
        else:
            obj.col_labels = col_labels;

        obj.weights = np.ones(data.shape);
        obj.filters = dict();

        return obj;
        
    def __init__(self, data, row_labels=[], col_labels=[]):
        pass;
                
    def __array_finalize__(self, obj):
        if obj is None:
            return;

        self.filters = getattr(obj, 'filters'   , []);
        self.weights    = getattr(obj, 'weights'   , []);
        self.row_labels = getattr(obj, 'row_labels', []);
        self.col_labels = getattr(obj, 'col_labels', []);
        
    def distance_matrix(self, filter_name=None):
        
        dist_vec = distance.pdist(self.projection_data(filter_name).T, metric='euclidean');
        return distance.squareform(dist_vec);

    def eval_signature(self, signature):
        """For a signature, calculate a score against each sample in the data
        
        Parameters
        ----------
        signature :   Signature,
            expression data to evaluate the signature against

        Returns
        -------
        out : 1D ndarray, length = Num_Samples 
            result of evaluating this signature on each sample in data
        
        """
        
        sig_vector = signature.sig_indices(self.row_labels);

        signature_norm = np.sum(np.abs(sig_vector));
        if(signature_norm == 0): #If no genes match signature
            raise ValueError("No genes match signature");

        weights = self.weights;
        pdata = self * sig_vector * weights;
        
        sig_scores = pdata.sum(axis=0);
        sig_scores /= np.sum(np.abs(sig_vector)*weights, axis=0); #Only normalize by weights in the signature

        sig_scores = scipy.stats.rankdata(sig_scores, method="average");
        
        return sig_scores;
        
    def subset_genes(self, indices):
        if(issubclass(type(indices), np.ndarray)):
            if(indices.dtype == np.bool):
                indices = np.nonzero(indices)[0];
        
        out = self[indices,:];
        out.weights = out.weights[indices,:];
        out.row_labels = [self.row_labels[i] for i in indices];

        return(out);
    
    def subset_samples(self, indices):
        if(issubclass(type(indices), np.ndarray)):
            if(indices.dtype == np.bool):
                indices = np.nonzero(indices)[0];
        
        out = self[:, indices];
        out.weights = out.weights[:, indices];
        out.col_labels = [self.col_labels[i] for i in indices];
        return(out);

    def projection_data(self, filter_name=None):
        """
        Returns the data matrix in self filtered by filters[filter_name]
        :return: A sliced copy of the underlying data matrix
        """
        if(filter_name is None or filter_name == 'None'):
            return self.base;
        else:
            filter = self.filters[filter_name];
            filter_i = [i for i,gene in enumerate(self.row_labels) if gene in filter];
            return self.base[filter_i,:];

    def projection_weights(self, filter_name=None):
        """
        Returns the weight matrix, filtered by projection_mask
        :return: A sliced copy of the weight matrix
        """

        if(filter_name is None or filter_name == 'None'):
            return self.weights;
        else:
            filter = self.filters[filter_name];
            filter_i = [i for i,gene in enumerate(self.row_labels) if gene in filter];
            return self.weights[filter_i,:];

class ProbabilityData(np.ndarray):
    
    def __new__(subtype, data, expression_data):
        
        obj = np.asarray(data).view(subtype);
        
        obj.row_labels = expression_data.row_labels;
        
        obj.col_labels = expression_data.col_labels;
        
        obj.expression_data = expression_data;

        obj.weights = expression_data.weights;

        return obj;
                
    def __array_finalize__(self, obj):
        if obj is None:
            return;

        self.weights    = getattr(obj, 'weights'   , []);
        self.row_labels = getattr(obj, 'row_labels', []);
        self.col_labels = getattr(obj, 'col_labels', []);
        self.expression_data = getattr(obj, 'expression_data', []);
        
    def distance_matrix(self, filter_name=None):

        pos_self = self.projection_data(filter_name)
        neg_self = 1-pos_self;
        
        prob_dist = np.dot(neg_self.T, pos_self) + \
                    np.dot(pos_self.T, neg_self);
        
        #Set all diagonal entries equal to zero
        np.fill_diagonal(prob_dist, 0.0);
        
        return prob_dist;

    #All data values in PCData have the same weight
    @property
    def filters(self):
        return self.expression_data.filters;

    def eval_signature(self, signature):
        """For a signature, calculate a score against each sample in the data
        
        Parameters
        ----------
        signature :   Signature,
            expression data to evaluate the signature against

        Returns
        -------
        out : 1D ndarray, length = Num_Samples 
            result of evaluating this signature on each sample in data
        
        """
        sig_vector = signature.sig_indices(self.row_labels);

        signature_norm = np.sum(np.abs(sig_vector));
        if(signature_norm == 0): #No genes match signature
            raise ValueError("No genes match signature");

        #Probability data already incorporates false-negative probability so no weights are used.
        weights = np.ones(self.shape);

        pdata = self * sig_vector * weights;
        
        sig_scores = pdata.sum(axis=0);
        sig_scores /= np.sum(np.abs(sig_vector)*weights, axis=0); #Only normalize by weights in the signature

        sig_scores = scipy.stats.rankdata(sig_scores, method="average");

        return sig_scores;
        
    def subset_genes(self, indices):
        if(issubclass(type(indices), np.ndarray)):
            if(indices.dtype == np.bool):
                indices = np.nonzero(indices)[0];
        
        out = self[indices,:];
        out.weights = out.weights[indices,:];
        out.row_labels = [self.row_labels[i] for i in indices];
        out.expression_data = out.expression_data.subset_genes(indices);
        return(out);
    
    def subset_samples(self, indices):
        if(issubclass(type(indices), np.ndarray)):
            if(indices.dtype == np.bool):
                indices = np.nonzero(indices)[0];
        
        out = self[:, indices];
        out.weights = out.weights[:, indices];
        out.col_labels = [self.col_labels[i] for i in indices];
        out.expression_data = out.expression_data.subset_samples(indices);
        return(out);

    def projection_data(self, filter_name=None):
        """
        Returns the data matrix in self filtered by filters[filter_name]
        :return: A sliced copy of the underlying data matrix
        """
        if(filter_name is None or filter_name == 'None'):
            return self.base;
        else:
            filter = self.filters[filter_name];
            filter_i = [i for i,gene in enumerate(self.row_labels) if gene in filter];
            return self.base[filter_i,:];

    def projection_weights(self, filter_name=None):
        """
        Returns the weight matrix, filtered by projection_mask
        :return: A sliced copy of the weight matrix
        """

        if(filter_name is None or filter_name == 'None'):
            return self.weights;
        else:
            filter = self.filters[filter_name];
            filter_i = [i for i,gene in enumerate(self.row_labels) if gene in filter];
            return self.weights[filter_i,:];

class PCData(np.ndarray):
    
    def __new__(subtype, data, variance, parent_data):
        
        obj = np.asarray(data).view(subtype);

        obj.variance = variance;
        
        obj.row_labels = ["PC"+str(i+1) for i in range(data.shape[0])];
        
        obj.col_labels = parent_data.col_labels;
        
        obj.parent_data = parent_data;
        
        return obj;

    #All data values in PCData have the same weight
    @property
    def weights(self):
        return np.ones(self.shape);

    def __array_finalize__(self, obj):
        if obj is None:
            return;

        self.variance = getattr(obj, 'variance', []);
        self.row_labels = getattr(obj, 'row_labels', []);
        self.col_labels = getattr(obj, 'col_labels', []);
        self.parent_data = getattr(obj, 'parent_data', []);
        
    def distance_matrix(self, filter_name=None):
        
        dist_vec = distance.pdist(self.T, metric='euclidean');
        return distance.squareform(dist_vec);

    def eval_signature(self, signature):
        """For a signature, calculate a score against each sample in the data
        
        Parameters
        ----------
        signature :   Signature,
            expression data to evaluate the signature against

        Returns
        -------
        out : 1D ndarray, length = Num_Samples 
            result of evaluating this signature on each sample in data
        
        """
        
        #On Principle Components, just evaluate signature on parent data
        return self.parent_data.eval_signature(signature);
        
    
    def subset_samples(self, indices):
        if(issubclass(type(indices), np.ndarray)):
            if(indices.dtype == np.bool):
                indices = np.nonzero(indices)[0];
        
        out = self[:, indices];
        out.col_labels = [self.col_labels[i] for i in indices];
        out.parent_data = out.parent_data.subset_samples(indices);
        return(out);

    def subset_components(self, indices):
        if(issubclass(type(indices), np.ndarray)):
            if(indices.dtype == np.bool):
                indices = np.nonzero(indices)[0];

        out = self[indices,:];
        out.variance = out.variance[indices];
        out.row_labels = [self.row_labels[i] for i in indices];
        return(out);


    def projection_data(self, filter_name=None):
        """
        Returns the data matrix in self filtered by projection_mask
        :return: A sliced copy of self
        """
        #PCA Data does not have projection_mask, just return all
        return self;

    def projection_weights(self, filter_name=None):
        """
        Returns the weight matrix, filtered by projection_mask
        :return: A sliced copy of the weight matrix
        """
        #PCA Data does not have a projection mask, just return weights
        return self.weights;
