# -*- coding: utf-8 -*-
"""
Created on Thu Mar 05 20:13:40 2015

@author: daved_000
"""

import numpy as np;
from scipy.spatial import distance;

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
        
        return obj;
        
    def __init__(self, data, row_labels=[], col_labels=[]):
        pass;
                
    def __array_finalize__(self, obj):
        if obj is None:
            return;
        
        self.row_labels = getattr(obj, 'row_labels', []);
        self.col_labels = getattr(obj, 'col_labels', []);
        
    def distance_matrix(self):
        
        dist_vec = distance.pdist(self.T, metric='euclidean');
        return distance.squareform(dist_vec);

    def eval_signature(self, signature, fn_prob=[]):
        """For a signature, calculate a score against each sample in the data
        
        Parameters
        ----------
        signature :   Signature,
            expression data to evaluate the signature against
        fn_prob : array-like, shape (Num_Genes, Num_Samples)
            false-negative probability used to weight signature score
            (Optional)
        
        Returns
        -------
        out : 1D ndarray, length = Num_Samples 
            result of evaluating this signature on each sample in data
        
        """
        
        sig_vector = signature.sig_indices(self.row_labels);

        signature_norm = np.sum(np.abs(sig_vector));
        if(signature_norm == 0): #If no genes match signature
            raise ValueError("No genes match signature");

        fn_prob = np.array(fn_prob);
        if(len(fn_prob)==0):
            weights = np.ones(self.shape);
        else:
            weights = 1-fn_prob;

        pdata = self * sig_vector * weights;
        
        sig_scores = pdata.sum(axis=0);
        sig_scores /= np.sum(np.abs(sig_vector)*weights, axis=0); #Only normalize by weights in the signature
        
        return sig_scores;
        
    def subset_genes(self, indices):
        if(type(indices) is np.ndarray):
            if(indices.dtype == np.bool):
                indices = np.nonzero(indices);
        
        out = self[indices,:];
        out.row_labels = [self.row_labels[i] for i in indices];
        return(out);
    
    def subset_samples(self, indices):
        if(type(indices) is np.ndarray):
            if(indices.dtype == np.bool):
                indices = np.nonzero(indices);
        
        out = self[:, indices];
        out.col_labels = [self.col_labels[i] for i in indices];
        return(out);
        
class ProbabilityData(np.ndarray):
    
    def __new__(subtype, data, expression_data):
        
        obj = np.asarray(data).view(subtype);
        
        obj.row_labels = expression_data.row_labels;
        
        obj.col_labels = expression_data.col_labels;
        
        obj.expression_data = expression_data;
        
        return obj;
                
    def __array_finalize__(self, obj):
        if obj is None:
            return;
        
        self.row_labels = getattr(obj, 'row_labels', []);
        self.col_labels = getattr(obj, 'col_labels', []);
        self.expression_data = getattr(obj, 'expression_data', []);
        
    def distance_matrix(self):
        
        neg_self = 1-self;
        
        prob_dist = np.dot(neg_self.T, self) + np.dot(self.T, neg_self);       
        
        #Set all diagonal entries equal to zero
        i = np.arange(prob_dist.shape[0]);
        prob_dist[i,i] = 0;        
        
        return prob_dist;

    def eval_signature(self, signature, fn_prob=[]):
        """For a signature, calculate a score against each sample in the data
        
        Parameters
        ----------
        signature :   Signature,
            expression data to evaluate the signature against
        fn_prob : array-like, shape (Num_Genes, Num_Samples)
            false-negative probability used to weight signature score
            (Optional). This is ignored for probability data as the probability
            values already capture the false-negative probability.
        
        Returns
        -------
        out : 1D ndarray, length = Num_Samples 
            result of evaluating this signature on each sample in data
        
        """
        sig_vector = signature.sig_indices(self.row_labels);

        signature_norm = np.sum(np.abs(sig_vector));
        if(signature_norm == 0): #No genes match signature
            raise ValueError("No genes match signature");

        weights = np.ones(self.shape);

        pdata = self * sig_vector * weights;
        
        sig_scores = pdata.sum(axis=0);
        sig_scores /= np.sum(np.abs(sig_vector)*weights, axis=0); #Only normalize by weights in the signature

        return sig_scores;
        
    def subset_genes(self, indices):
        if(type(indices) is np.ndarray):
            if(indices.dtype == np.bool):
                indices = np.nonzero(indices);
        
        out = self[indices,:];
        out.row_labels = [self.row_labels[i] for i in indices];
        out.expression_data = out.expression_data.subset_genes(indices);
        return(out);
    
    def subset_samples(self, indices):
        if(type(indices) is np.ndarray):
            if(indices.dtype == np.bool):
                indices = np.nonzero(indices);
        
        out = self[:, indices];
        out.col_labels = [self.col_labels[i] for i in indices];
        out.expression_data = out.expression_data.subset_samples(indices);
        return(out);
        
class PCData(np.ndarray):
    
    def __new__(subtype, data, parent_data):
        
        obj = np.asarray(data).view(subtype);
        
        obj.row_labels = ["PC"+str(i+1) for i in range(data.shape[0])];
        
        obj.col_labels = parent_data.col_labels;
        
        obj.parent_data = parent_data;
        
        return obj;
                
    def __array_finalize__(self, obj):
        if obj is None:
            return;
        
        self.row_labels = getattr(obj, 'row_labels', []);
        self.col_labels = getattr(obj, 'col_labels', []);
        self.parent_data = getattr(obj, 'parent_data', []);
        
    def distance_matrix(self):
        
        dist_vec = distance.pdist(self.T, metric='euclidean');
        return distance.squareform(dist_vec);

    def eval_signature(self, signature, fn_prob=[]):
        """For a signature, calculate a score against each sample in the data
        
        Parameters
        ----------
        signature :   Signature,
            expression data to evaluate the signature against
        fn_prob : array-like, shape (Num_Genes, Num_Samples)
            false-negative probability used to weight signature score
            (Optional)
        
        Returns
        -------
        out : 1D ndarray, length = Num_Samples 
            result of evaluating this signature on each sample in data
        
        """
        
        #On Principle Components, just evaluate signature on parent data
        return self.parent_data.eval_signature(signature, fn_prob);
        
    
    def subset_samples(self, indices):
        if(type(indices) is np.ndarray):
            if(indices.dtype == np.bool):
                indices = np.nonzero(indices);
        
        out = self[:, indices];
        out.col_labels = [self.col_labels[i] for i in indices];
        out.parent_data = out.parent_data.subset_samples(indices);
        return(out);

    def subset_components(self, indices):
        if(type(indices) is np.ndarray):
            if(indices.dtype == np.bool):
                indices = np.nonzero(indices);

        out = self[indices,:];
        out.row_labels = [self.row_labels[i] for i in indices];
        return(out);
