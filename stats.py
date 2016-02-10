import numpy as np

def print_dict_stats_with_ref(data,ref,keys,title=""):
    print "\n"+ title
    print "\n{:>12}  {:<10} {:<10} {:<10} {:<10} {:<10}".format(
                                                "","ref","mean","std","95 decile","95-mean")
    for i,key in enumerate(keys):
        print "{:>12} {:>10.3e} {:>10.3e} {:>10.3e} {:>10.3e} {:>10.3e}".format(key,
                        ref[key],
                        data[key].mean(),
                        data[key].std(),
                        np.percentile(data[key],q=95),
                        np.percentile(data[key],q=95) - data[key].mean())
