import numpy as np
import pyqt_fit as fitter
import pyqt_fit.nonparam_regression as smooth

def print_dict_stats_with_ref(data,ref,keys,title="",q=95):
    print("\n"+ title)
    print("\n{:>12}  {:<10} {:<10} {:<10} {:<10}".format(
                                                "","ref","mean","std/mean","{} decile".format(q)))
    for i,key in enumerate(keys):
        print("{:>12} {:>10.5f} {:>10.4f} {:>10.4f} {:>10.4f}".format(key,
                        ref[key],
                        data[key].mean(),
                        data[key].std()/data[key].mean(),
                        np.percentile(data[key],q=q)))

def rmse_at_zero_of_topq(data,q,ax=0):
    zero = np.zeros(data.shape)
    data = np.where(data > np.percentile(data,q=q),data,zero)
    return np.sqrt( np.mean(data**2,axis=ax) )

def kern_reg(xs,ys,cutx=8,bins=63):
    reg = smooth.NonParamRegression(xs,ys,
                                method=fitter.npr_methods.SpatialAverage())
    reg.fit()
    count, evxs = np.histogram(xs,bins)
    tmp, = np.where(count > cutx)
    evxs = evxs[slice(tmp[0],tmp[-1])]
    return evxs,reg.evaluate(evxs)

