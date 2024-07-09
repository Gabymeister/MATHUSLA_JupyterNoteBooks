import os
import numpy as np
## The Python Standard Library ##
import os,sys
import re
import ast
import copy as cp
from glob import glob
import types
from tqdm import tqdm
from types import ModuleType, FunctionType
from gc import get_referents
from array import array
import warnings
import operator
import inspect
import itertools
import functools
from functools import reduce # import needed for python3; builtin in python2
from collections import defaultdict
from datetime import datetime # Pylab will import datatime directly, so need to put this line after pylab..



import pickle
import joblib
import numpy as np
import pandas as pd
import scipy
import scipy.signal as signal
from scipy.io import loadmat
from scipy import stats,signal
from scipy import optimize
from scipy.optimize import curve_fit
from sklearn.cluster import KMeans
import uncertainties
import h5py

import ROOT as root
from root_numpy import array2tree

def most_frequent(nums):
    
    """Get the most frequent element in a list"""
    return max(set(nums), key = nums.count) 

def unzip(concatlist, divider=-1):
    lists = []
    n = 0
    for val in concatlist:
        if val == -1:
            n += 1
        else:
            while len(lists) <= n:
                lists.append([])
            lists[n].append(val)
    return lists

# Verification:
# print(coord_det2cms([0,0,1]), coord_cms2det(coord_det2cms([0,0,1])))
# tracks=get_truthtrack(ev)

def theta2eta(theta):
    return -np.log(np.tan(theta/2))


def quadsum(data):
    return np.sqrt(np.sum(np.power(data,2)))

def c2list(cvector):
    return [cvector.at(i) for i in range(cvector.size())]

def sortbyrow(a, irow):
    return a[:,a[irow, :].argsort()]
def sortbycol(a, icol):
    return a[a[:, icol].argsort()]


class Utils:
    class color:
        PURPLE = '\033[95m'
        CYAN = '\033[96m'
        DARKCYAN = '\033[36m'
        BLUE = '\033[94m'
        GREEN = '\033[92m'
        YELLOW = '\033[93m'
        RED = '\033[91m'
        BOLD = '\033[1m'
        UNDERLINE = '\033[4m'
        END = '\033[0m'
        
        @classmethod
        def purple(cls,ustring):
            return cls.PURPLE + ustring + cls.END
        @classmethod
        def cyan(cls,ustring):
            return cls.CYAN + ustring + cls.END
        @classmethod
        def darkcyan(cls,ustring):
            return cls.DARKCYAN + ustring + cls.END
        @classmethod
        def blue(cls,ustring):
            return cls.BLUE + ustring + cls.END
        @classmethod
        def green(cls,ustring):
            return cls.GREEN + ustring + cls.END
        @classmethod
        def yellow(cls,ustring):
            return cls.YELLOW + ustring + cls.END
        @classmethod
        def red(cls,ustring):
            return cls.RED + ustring + cls.END
        @classmethod
        def bold(cls,ustring):
            return cls.BOLD + ustring + cls.END
        @classmethod
        def underline(cls,ustring):
            return cls.UNDERLINE + ustring + cls.END

    @staticmethod
    def groupby(key,item):
        """
        Group the second list `item` with the first list `key`

        Returns
        -------
        key_grouped
        item_grouped
        """
        data = [(item_i,key_i) for key_i, item_i in zip(key,item)]
        res = defaultdict(list)
        for v, k in data: res[k].append(v)

        key_grouped = [k for k,v in res.items()]
        item_grouped = [v for k,v in res.items()]

        return key_grouped,item_grouped


    @staticmethod
    def is_odd(num):
        return num & 0x1

    @staticmethod
    def getsize(obj):
        """sum size of object & members."""
        # Custom objects know their class.
        # Function objects seem to know way too much, including modules.
        # Exclude modules as well.

        BLACKLIST = type, ModuleType, FunctionType
        if isinstance(obj, BLACKLIST):
            raise TypeError('getsize() does not take argument of type: '+ str(type(obj)))
        seen_ids = set()
        size = 0
        objects = [obj]
        while objects:
            need_referents = []
            for obj in objects:
                if not isinstance(obj, BLACKLIST) and id(obj) not in seen_ids:
                    seen_ids.add(id(obj))
                    size += sys.getsizeof(obj)
                    need_referents.append(obj)
            objects = get_referents(*need_referents)
        return size

    @staticmethod
    def find_float(string,return_format="float"):
        """
        Find all float numbers inside a string

        return_format: "float" or "str"
        """
        float_in_str = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", string)
        if return_format == "float":
            data = [float(i) for i in float_in_str]
        elif return_format=="str":
            data = [i for i in float_in_str]
        return data

    @staticmethod
    def strike(text):
        return ''.join([u'\u0336{}'.format(c) for c in text])

    @staticmethod
    def getIVserlist(ser0,nser,rawdir='/gpfs/slac/staas/fs1/g/supercdms/tf/nexus/midasdata/NEXUS/R7/Raw/'):
        """
        Get a list of series after an initial series
        """
        fnames = sorted([fname for fname in os.listdir(rawdir) ])
        return fnames[fnames.index(ser0):fnames.index(ser0)+nser]

    @staticmethod
    def append_dicts(dict1,dict2,verbose=False):
        """
        Append dict2 after dict1
        """
        dict_combined=dict()
        for key in dict1:
            #dict_combined[key]=np.append(dict1[key],dict2[key])
            if key in dict2:
                try:
                    dict_combined[key]=np.concatenate((dict1[key],dict2[key]))#.astype(type(dict2[key][0]))
                except Exception as e:
                    if verbose:
                        print(key,e)
                    continue
        return dict_combined

    @staticmethod
    def flatten1d(a):
        #functools_reduce_iconcat
        """Reduce an array to 1-d"""
        return np.array(functools.reduce(operator.iconcat, a, []))

    @staticmethod
    def center(x):
        return 0.5*(x[1:]+x[:-1])
    
    @staticmethod
    def contrast(a,b,weight=1):
        """return (a-b*weight)/(a+b*weight)
        """
        return (a-b*weight)/(a+b*weight)

    @staticmethod
    def pch(inner,outer,scale=1):
        """return (outer*(1+scale)-inner*(1-scale))/(outer*(1+scale)+inner*(1-scale))
        """
        return (outer*(1+scale)-inner*(1-scale))/(outer*(1+scale)+inner*(1-scale))


    @staticmethod
    def parabola_3points(x1,x2,x3,y1,y2,y3,return_vertex=False):
        """
        Solve for A,B,C for parabola: A x1^2 + B x1 + C = y1 with 3 points
        """
        denom = (x1-x2) * (x1-x3) * (x2-x3);
        A     = (x3 * (y2-y1) + x2 * (y1-y3) + x1 * (y3-y2)) / denom;
        B     = (x3*x3 * (y1-y2) + x2*x2 * (y3-y1) + x1*x1 * (y2-y3)) / denom;
        C     = (x2 * x3 * (x2-x3) * y1+x3 * x1 * (x3-x1) * y2+x1 * x2 * (x1-x2) * y3) / denom;

        x_vertex = -B / (2*A);
        y_vertex = C - B*B / (4*A);
        if return_vertex:
            return x_vertex,y_vertex

        return A,B,C



    @staticmethod
    def calc_parabola_vertex(x1, y1, x2, y2, x3, y3):
        '''
        3-point parabola:
        http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
        '''

        denom = (x1-x2) * (x1-x3) * (x2-x3);
        A     = (x3 * (y2-y1) + x2 * (y1-y3) + x1 * (y3-y2)) / denom;
        B     = (x3*x3 * (y1-y2) + x2*x2 * (y3-y1) + x1*x1 * (y2-y3)) / denom;
        C     = (x2 * x3 * (x2-x3) * y1+x3 * x1 * (x3-x1) * y2+x1 * x2 * (x1-x2) * y3) / denom;

        return A,B,C

    @staticmethod
    def parabola_interpolation(data,ind):
        """
        parabola interpolation to find the y value at a given x, with 3 points in data around a index (x)

        Inputs
        -------
        data: list of y value
        ind: x value
        """
        if type(ind) is int:
            return data[ind]
        else:
            xa=int(np.floor(ind)); xa=min(xa,len(data)-2)
            xb=int(np.ceil(ind));  xb=min(xb,len(data)-1)
            ya,yb=(data[xa],data[xb])         
            if ind > (xa+xb)/2:
                if xb+1 < len(data):
                    xc = xb+1
                else:
                    xc = xa-1
            else:
                if xa-1 >= 0:
                    xc = xa-1
                else:
                    xc = xb+1
            yc = data[xc]
            A,B,C = Utils.parabola_3points(xa, xb, xc, ya, yb, yc)
            y_vertex = C - B*B / (4*A)
            return y_vertex

    @staticmethod
    def linear_interpolation(data,ind):
        """
        Linear interpolation to find the y value at a given x, with 2 points in data around a index (x)

        Inputs
        -------
        data: list of y value
        ind: x value
        """
        if type(ind) is int:
            return data[ind]
        else:
            xa=int(np.floor(ind)); xa=min(xa,len(data)-2)
            xb=int(np.ceil(ind));  xb=min(xb,len(data)-1)
            ya,yb=(data[xa],data[xb])
            xc=ind
            m = (ya - yb) / (xa - xb)
            yc = (xc - xb) * m + yb
            return yc

    @staticmethod
    def bilinear_interpolation(x, y, points):
        '''Interpolate (x,y) from values associated with four points.

        The four points are a list of four triplets:  (x, y, value).
        The four points can be in any order.  They should form a rectangle.

            >>> bilinear_interpolation(12, 5.5,
            ...                        [(10, 4, 100),
            ...                         (20, 4, 200),
            ...                         (10, 6, 150),
            ...                         (20, 6, 300)])
            165.0

        '''
        # See formula at:  http://en.wikipedia.org/wiki/Bilinear_interpolation

        points = sorted(points)               # order points by x, then by y
        (x1, y1, q11), (_x1, y2, q12), (x2, _y1, q21), (_x2, _y2, q22) = points

        if x1 != _x1 or x2 != _x2 or y1 != _y1 or y2 != _y2:
            raise ValueError('points do not form a rectangle')
        if not x1 <= x <= x2 or not y1 <= y <= y2:
            raise ValueError('(x, y) not within the rectangle')

        return (q11 * (x2 - x) * (y2 - y) +
                q21 * (x - x1) * (y2 - y) +
                q12 * (x2 - x) * (y - y1) +
                q22 * (x - x1) * (y - y1)
               ) / ((x2 - x1) * (y2 - y1) + 0.0)

    @staticmethod
    def bilinear_interpolation_array(x, y, data):
        if type(x) is int and type(y) is int:
            return(data[x,y])
        x_l,x_h = int(np.floor(x)),int(np.ceil(x))
        y_l,y_h   = int(np.floor(y)),int(np.ceil(y))
        if (x_l==x_h):
            return Utils.linear_interpolation(data[x_l,:],y)
        if (y_l==y_h):
            return Utils.linear_interpolation(data[:,y_l],x)
        return Utils.bilinear_interpolation(x, y, [(x_l,y_l,data[x_l,y_l]),(x_l,y_h,data[x_l,y_h%data.shape[1]]),(x_h,y_l,data[x_h%data.shape[0],y_l]),(x_h,y_h,data[x_h%data.shape[0],y_h%data.shape[1]])])

    @staticmethod
    def lin_interp(x, y, i, half):
        """
        Linear interpolation to find the x value at a given y, with 2 points in data around a index (i) and the y value (half)

        Inputs
        -------
        data: list of y value
        ind: x value
        """
        return x[i] + (x[i+1] - x[i]) * ((half - y[i]) / (y[i+1] - y[i]))

    @staticmethod
    def fwhm(x, y, height=0.5):
        """
        Find the interval in x corrsponding to full width half maximum
        x: list of x values
        y: list of y values
        height: default to 0.5 (half maximum)
        """
        half = max(y)*height
        signs = np.sign(np.add(y, -half))
        zero_crossings = (signs[0:-2] != signs[1:-1])
        zero_crossings_i = np.where(zero_crossings)[0]
        return [Utils.lin_interp(x, y, zero_crossings_i[0], half),
                Utils.lin_interp(x, y, zero_crossings_i[1], half)]

    @staticmethod
    def find_crossing(x, y, height=1):
        """
        Find the interval in x corrsponding to a certain value in y
        x: list of x values
        y: list of y values
        height: default to 1
        """
        half = height
        signs = np.sign(np.add(y, -half))
        zero_crossings = (signs[0:-2] != signs[1:-1])
        zero_crossings_i = np.where(zero_crossings)[0]
        return [Utils.lin_interp(x, y, zero_crossings_i[0], half),
                Utils.lin_interp(x, y, zero_crossings_i[1], half)]

    #----------------------------------------------------------------------------
    # A few math functions redefined

    @staticmethod
    def nll_poisson(data, model):
        nll = np.sum([-scipy.stats.poisson.logpmf(k,mu) for k,mu in zip(data,model)])
        return nll

    @staticmethod
    def Uniform(x,A):
        if type(x) not in [np.ndarray,list]:
            return A
        else:
            return A*np.ones_like(x)

    @staticmethod
    def Exp(x,A,t):
        return A*np.exp(-x/t)

    @staticmethod
    def Gauss(x, A, mean, sigma):
        return A * np.exp(-(x - mean)**2 / (2 * sigma**2))

    @staticmethod
    def Gauss_sideband(x, A, mean, sigma, a1,a2):
        # a1 for left, a2 for right
        return Utils.Gauss(x, A, mean, sigma) + sqrt(2*np.pi)*sigma/2*(a1*scipy.special.erfc((x-mean)/sqrt(2)/sigma) + a2*(2-scipy.special.erfc((x-mean)/sqrt(2)/sigma)))

    @staticmethod
    def Poisson(k, Lambda,A):
        # Lambda: mean, A: amplitude
        return A*(Lambda**k/scipy.special.factorial(k)) * np.exp(-Lambda)

    @staticmethod
    def Poly(x, *P):
        '''
        Compute polynomial P(x) where P is a vector of coefficients
        Lowest order coefficient at P[0].  Uses Horner's Method.
        '''
        result = 0
        for coeff in P[::-1]:
            result = x * result + coeff
        return result

    @staticmethod
    def Chi2_reduced(x, dof, A):
        return scipy.stats.chi2.pdf(x*dof,dof)*A

    
    @staticmethod
    def fitu(array, fit_range=None, n_bins=1000, functions=(root.RooGaussian,),
            initial_values=((0.0, 1.0, 1.0),),
            bounds=(((-1e6, 1e6), (0, 1e6), (0, 1e6)), ),
            set_constant=None,
            verbosity=0):
        """Uses the RooFit package to fit a dataset (instead of fitting a histogram)
        Source: Sasha Zaytsev

        Parameters
        ----------
        array : 1-d array or list
            input data array to fit
        fit_range : tuple
            data range for the fit (x_lower, x_upper)
        n_bins : int
            number of points on the x-axis in the output. Does not affect the fit!
        functions : tuple of RooAbsPdf
            Roo pdf function.
            Examples:
            RooGaussian, RooUniform, RooPolynomial, RooExponential
        initial_values : tuple of tuples of floats
            inital values of parameters
            Example:
            functions=(root.RooGaussian, root.RooExponential, root.Polynomial), initial_values=((mean, sigma, a), (exp_k, exp_a), (p1, p2, ..., a))
        bounds : tuple of tuples of tuples of floats
            min and max allowed parameter values
            Example:
            functions=(root.RooGaussian, root.RooExponential), bounds=(((min_mean, max_mean),(min_sig,max_sig),(min_a, max_a)), ((min_k, max_k),(min_a, max_a)))
        set_constant : tuple of tuples of bools   or   None
            whether to fix a certain parameter at a constant value.
            If equals to None, then none of the parameters is fixed
            Example:
            functions=(root.RooGaussian, root.RooExponential), set_constant=((fix_mean, fix_sigma), (fix_k))
        verbosity : int
            verbosity level (might not work. It's tricky)
            -2 - print nothing
            -1 - print errors
             0 - print errors and fit results
             1 - print warnings
             2 - print info

        Returns
        -------
        x, y, param_values, param_errors
        x : array
            bin centers
        y : array
            fit function values
        param_values : tuple of tuples
            values of fitted parameters. Has the same shape as 'initial_values' arg
        param_values : tuple of tuples
            errors of fitted parameters. Has the same shape as 'initial_values' arg
        """

        # trying to suppress output
        if verbosity < -1:
            root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.FATAL)
        if verbosity == -1 or verbosity == 0:
            root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
        if verbosity == 1:
            root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
        if verbosity >= 2:
            root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.INFO)

        if type(array)==list:
            array = np.array(array)

        if fit_range is None:
            fit_range = (np.min(array), np.max(array))

        # create a tree with one branch
        tree = array2tree(np.array(array, dtype=[('data', np.float64)]))

        data_var = root.RooRealVar('data', 'data', fit_range[0], fit_range[1])
        data_arg_set = root.RooArgSet(data_var)

        dataset = root.RooDataSet('dataset', 'dataset', tree, data_arg_set)

        parameters = []
        roo_functions = []
        amplitudes = []

        # iterating through the functions
        func_names = []
        for i,f in enumerate(functions):
            func_name = f.__name__
            # remove the Roo prefix
            if len(func_name)>3 and func_name[:3]=='Roo':
                func_name = func_name[3:]

            base_func_name = func_name
            k = 2
            while func_name in func_names:
                func_name = '%s%i'%(base_func_name, k)
                k+=1

            func_names.append(func_name)

            # creating function parameters
            func_parameters = []
            for j,initial_value in enumerate(initial_values[i][:-1]):
                name = '%s_p%i'%(func_name, j)
                parameter = root.RooRealVar(name, name, initial_value, *bounds[i][j])
                if not(set_constant is None) and set_constant[i][j]:
                    parameter.setConstant(True)
                func_parameters.append(parameter)
            parameters.append(func_parameters)

            # creating function amplitude
            name = '%s_a'%(func_name)
            amplitudes.append(root.RooRealVar(name, name, initial_values[i][-1], *bounds[i][-1]))

            if func_name=='Polynomial':
                roo_functions.append(f(func_name, func_name, data_var, root.RooArgList(*func_parameters)))
            elif func_name=='Uniform' or len(func_parameters)==0:
                roo_functions.append(f(func_name, func_name, data_arg_set))
            else:
                roo_functions.append(f(func_name, func_name, data_var, *func_parameters))

        function_list = root.RooArgList(*roo_functions)
        amplitude_list = root.RooArgList(*amplitudes)
        pdf = root.RooAddPdf('pdf', 'pdf', function_list, amplitude_list)

        # fitting
        fit_results = pdf.fitTo(dataset, root.RooFit.Save(), root.RooFit.Range(*fit_range), root.RooFit.PrintLevel(verbosity-1))
        if fit_results.status()!=0:
            if verbosity>=-1:
                print('----- FIT STATUS != 0 -----')
        if verbosity>=0:
            fit_results.Print()

        tf_parameters = []
        param_values = []
        param_errors = []

        for i,params in enumerate(parameters):
            tf_parameters += params
            param_values.append([p.getVal() for p in params] + [amplitudes[i].getVal()])
            param_errors.append([p.getError() for p in params] + [amplitudes[i].getError()])

        tf_parameters += amplitudes

        tf = pdf.asTF(root.RooArgList(data_var), root.RooArgList(*tf_parameters), data_arg_set)
        a = 0
        for amplitude in amplitudes:
            a += amplitude.getVal()

        bin_w = (fit_range[1] - fit_range[0])/n_bins
        x = np.linspace(fit_range[0]+bin_w/2, fit_range[1]-bin_w/2, n_bins)
        y = np.array([a*tf.Eval(ix) for ix in x])*bin_w

        return x, y, param_values, param_errors
    
    
def pull(x_measure, x_truth, x_unc):
    return (x_measure-x_truth)/x_unc

def poissonerror_div(N1,N2):
    return np.sqrt(1/N1+1/N2)*N1/N2

def chi2_calc(x_est, x_true, err):
    return sum([(x_est[i]-x_true[i])**2/err[i]**2 for i in range(len(err))])


