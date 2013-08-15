## Copyright (C) 2013 Lars Simon Zehnder
#
# This file is part of finmix.
#
# finmix is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# finmix is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Rcpp.  If not, see <http://www.gnu.org/licenses/>.

## Class 'model' --------------------------------------------------
setGeneric("mixturemar", function(object, J) standardGeneric("mixturemar"))

setGeneric("getDist", function(object) standardGeneric("getDist"))

setGeneric("getR", function(object) standardGeneric("getR"))

setGeneric("getK", function(object) standardGeneric("getK"))

setGeneric("getWeight", function(object) standardGeneric("getWeight"))

setGeneric("getPar", function(object) standardGeneric("getPar"))

setGeneric("getIndicmod", function(object) standardGeneric("getIndicmod"))

setGeneric("getIndicfix", function(object) standardGeneric("getIndicfix"))

setGeneric("getT", function(object) standardGeneric("getT"))

setGeneric("setDist<-", function(object, value) standardGeneric("setDist<-"))

setGeneric("setR<-", function(object, value) standardGeneric("setR<-"))

setGeneric("setK<-", function(object, value) standardGeneric("setK<-"))

setGeneric("setWeight<-", function(object, value) standardGeneric("setWeight<-"))

setGeneric("setPar<-", function(object, value) standardGeneric("setPar<-"))

setGeneric("setIndicmod<-", function(object, value) standardGeneric("setIndicmod<-"))

setGeneric("setIndicfix<-", function(object, value) standardGeneric("setIndicfix<-"))

setGeneric("setT<-", function(object, value) standardGeneric("setT<-"))

## Class 'modelmoments' --------------------------------------------

setGeneric("getMean", function(object) standardGeneric("getMean"))

setGeneric("getVar", function(object) standardGeneric("getVar"))

setGeneric("getModel", function(object) standardGeneric("getModel"))

## Class 'cmodelmoments' -------------------------------------------

setGeneric("getHigher", function(object) standardGeneric("getHigher"))

setGeneric("getSkewness", function(object) standardGeneric("getSkewness"))

setGeneric("getKurtosis", function(object) standardGeneric("getKurtosis"))

## Class 'dmodelmoments' -------------------------------------------

setGeneric("getOver", function(object) standardGeneric("getOver"))

setGeneric("getFactorial", function(object) standardGeneric("getFactorial"))

setGeneric("getZero", function(object) standardGeneric("getZero"))

## Class 'normultmodelmoments' -------------------------------------

setGeneric("generateMoments", function(object) standardGeneric("generateMoments"))

setGeneric("getB", function(object) standardGeneric("getB"))

setGeneric("getW", function(object) standardGeneric("getW"))

setGeneric("getRdet", function(object) standardGeneric("getRdet"))

setGeneric("getRtr", function(object) standardGeneric("getRtr"))

setGeneric("getCorr", function(object) standardGeneric("getCorr"))

## Class 'exponentialmodelmoments' ---------------------------------

setGeneric("getExtrabinvar", function(object) standardGeneric("getExtrabinvar"))

## Class 'data' ----------------------------------------------------

setGeneric("getY", function(object) standardGeneric("getY"))

setGeneric("getN", function(object) standardGeneric("getN"))

setGeneric("getS", function(object) standardGeneric("getS"))

setGeneric("getBycolumn", function(object) standardGeneric("getBycolumn"))

setGeneric("getName", function(object) standardGeneric("getName"))

setGeneric("getType", function(object) standardGeneric("getType"))

setGeneric("getSim", function(object) standardGeneric("getSim"))

setGeneric("getExp", function(object) standardGeneric("getExp"))

setGeneric("setY<-", function(object, value) standardGeneric("setY<-"))

setGeneric("setN<-", function(object,  value) standardGeneric("setN<-"))

setGeneric("setS<-", function(object, value) standardGeneric("setS<-"))

setGeneric("setBycolumn<-", function(object, value) standardGeneric("setBycolumn<-"))

setGeneric("setName<-", function(object, value) standardGeneric("setName<-"))

setGeneric("setType<-", function(object, value) standardGeneric("setType<-"))

setGeneric("setSim<-", function(object, value) standardGeneric("setSim<-"))

setGeneric("setExp<-", function(object, value) standardGeneric("setExp<-"))

## Class 'groupmoments' ----------------------------------------------

setGeneric("getNK", function(object) standardGeneric("getNK"))

setGeneric("getWK", function(object) standardGeneric("getWK"))

setGeneric("getData", function(object) standardGeneric("getData"))

## Class 'sdatamoments' ----------------------------------------------

setGeneric("getGmoments", function(object) standardGeneric("getGmoments"))

## Class 'csdatamoments' ---------------------------------------------

setGeneric("getSmoments", function(object) standardGeneric("getSmoments"))

## Class 'prior' -----------------------------------------------------

setGeneric("generatePrior", function(object, data, model, varargin, coef.mat)
           {
               standardGeneric("generatePrior")
           }
)

setGeneric("getHier", function(object) standardGeneric("getHier"))

setGeneric("setHier<-", function(object, value) standardGeneric("setHier<-"))

## Class 'mcmc' -------------------------------------------------------

setGeneric("getBurnin", function(object) standardGeneric("getBurnin"))

setGeneric("getM", function(object) standardGeneric("getM"))

setGeneric("getStartpar", function(object) standardGeneric("getStartpar"))

setGeneric("getStoreS", function(object) standardGeneric("getStoreS"))

setGeneric("getStorepost", function(object) standardGeneric("getStorepost"))

setGeneric("getRanperm", function(object) standardGeneric("getRanperm"))

setGeneric("setBurnin<-", function(object, value) standardGeneric("setBurnin<-"))

setGeneric("setM<-", function(object, value) standardGeneric("setM<-"))

setGeneric("setStartpar<-", function(object, value) standardGeneric("setStartpar<-"))

setGeneric("setStoreS<-", function(object, value) standardGeneric("setStoreS<-"))

setGeneric("setStorepost<-", function(object, value) standardGeneric("setStorepost<-"))

setGeneric("setRanperm<-", function(object, value) standardGeneric("setRanperm<-"))

















