# -*- coding: utf-8 -*-
"""Wrapper classes for different types of data

Here are the wrapper classes for ExpressionData, and PCData

This was created to organize nuances in how signature
scores, distance matrices, and anything else, is computed
on the different types of data.

"""

from __future__ import absolute_import, print_function, division;
import numpy as np;


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

    def filtered_genes(self, filter_name=None):
        """
        Returns the list of genes (in the order used by projection_data or projection_weights)
        filtered by the listed filter

        :param filter_name: Name of the filter to use
        :return: List of gene (row) names
        """

        if(filter_name is None or filter_name == 'None'):
            return self.row_labels;
        else:
            filter = self.filters[filter_name];
            fg = [gene for gene in self.row_labels if gene in filter]
            return fg;

    def projection_data(self, filter_name=None):
        """
        Returns the data matrix in self filtered by filters[filter_name]
        :return: A sliced copy of the underlying data matrix
        """
        if(filter_name is None or filter_name == 'None'):
            return self.base;
        else:
            filter = self.filters[filter_name];
            filter_i = [i for i, gene in enumerate(self.row_labels) if gene in filter];
            return self.base[filter_i, :];

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

    def get_normalized_copy(self, norm_method):
        """
        Runs norm_method on the data and returns an ExpressionData copy

        norm_method is of the form:

            out = norm_method(input);

            input: numpy.ndarray N x M
            out: numpy.ndarray N x M

        """

        normalized_data = norm_method(self.base);
        normalized_edata = ExpressionData(normalized_data,
                row_labels=self.row_labels,
                col_labels=self.col_labels);
        
        normalized_edata.weights = self.weights;
        normalized_edata.filters = self.filters;

        return normalized_edata;

    def merge_data(self, other):
        """
        Merges columns of self and other
        :param other: Instance of Expression Data with same rows as self
        :return: New Expression Data Instance combining self and other
        """

        #Ensure row labels match
        if(len(self.row_labels) != len(other.row_labels)):
            raise ValueError("Cant merge ExpressionData objects with different row labels");

        for x,y in zip(self.row_labels, other.row_labels):
            if(x != y):
                raise ValueError("Cant merge ExpressionData objects with different row labels");

        row_labels_merge = self.row_labels; #Could be either, they are ensured to be identical at this point

        data_merge = np.concatenate((self.base, other.base), axis=1);

        #For Filters, combine sets preferring self over other
        filters_merge = self.filters.copy();
        for key in other.filters:
            if(key not in filters_merge):
                filters_merge.update({key: other.filters[key]});

        weights_merge = np.concatenate((self.weights, other.weights), axis=1);
        col_labels_merge = self.col_labels + other.col_labels;

        merge_data = ExpressionData(data_merge, row_labels_merge, col_labels_merge);
        merge_data.weights = weights_merge;
        merge_data.filters = filters_merge;
        return merge_data;


class PCData(np.ndarray):
    
    def __new__(subtype, data, variance, loadings, parent_data):
        
        obj = np.asarray(data).view(subtype);

        obj.variance = variance;
        
        obj.row_labels = ["PC"+str(i+1) for i in range(data.shape[0])];

        obj.loadings = loadings;

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
        self.loadings = getattr(obj, 'loadings', []);
        self.parent_data = getattr(obj, 'parent_data', []);

#     def eval_signature(self, signature):
#         """For a signature, calculate a score against each sample in the data
#         
#         Parameters
#         ----------
#         signature :   Signature,
#             expression data to evaluate the signature against
# 
#         Returns
#         -------
#         out : 1D ndarray, length = Num_Samples 
#             result of evaluating this signature on each sample in data
#         
#         """
#         
#         #On Principle Components, just evaluate signature on parent data
#         return self.parent_data.eval_signature(signature);
        
    
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
        out.loadings = out.loadings[:,indices];
        out.variance = out.variance[indices];
        out.row_labels = [self.row_labels[i] for i in indices];
        return(out);

    def filtered_genes(self, filter_name=None):
        """
        Returns the list of genes (in the order used by projection_data or projection_weights)
        filtered by the listed filter

        :param filter_name: Name of the filter to use
        :return: List of gene (row) names
        """

        return self.parent_data.filtered_genes(filter_name);

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


