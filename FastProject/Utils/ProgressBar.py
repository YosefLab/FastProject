# -*- coding: utf-8 -*-
"""Simple progress bar used throughout FastProject

"""
from __future__ import absolute_import, print_function, division;

import sys

class ProgressBar:
    
    width = 40;  #Number of spaces to use    
    
    def __init__(self, max_i):
        self.max_i = max_i;
        self.i = 0;

        if(self.max_i == 0):
            self.max_i = 1;
        
        sys.stdout.write(self.build_str());
        sys.stdout.flush();
        
    def build_str(self):
        #Determine number of hashes
        n_hash = int(self.i / self.max_i * ProgressBar.width);
        
        #Determine percentage to display
        pcent = int(self.i / self.max_i * 100);
        pcent_str = '{:3d}'.format(pcent)
        
        out_string = '|';
        out_string = out_string + ''.join(['#']*n_hash);
        out_string = out_string + ''.join(['-']*(ProgressBar.width-n_hash));
        out_string = out_string + '| ' + pcent_str + '%'; 
        
        return out_string;
        
    def update(self):
        if(self.i < self.max_i):
            self.i = self.i + 1;

            #Only update stdout if number of hashes has changed
            if(int(self.i / self.max_i * ProgressBar.width) > int((self.i - 1) / self.max_i * ProgressBar.width)):

                new_str = self.build_str();

                sys.stdout.write(''.join(['\b']*len(new_str)));  #Erase old
                sys.stdout.write(new_str);

            if(self.i == self.max_i):
                sys.stdout.write('\n');
        
        sys.stdout.flush();

    def complete(self):
        #Only does anything if we haven't completed yet
        if(self.i != self.max_i):
            self.i = self.max_i-1;
            self.update();
