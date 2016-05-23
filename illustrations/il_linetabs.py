import pickle as pic

tabhead="""
\\begin{table}[]
\centering
    \\begin{tabular}{lccccl}
Species & $\lambda$  & E${}_l$ & $\log gf$ & depth & Comment\\\\ \hline"""
tabline = "{:<10}  & {:7.4f}   & {} & {}  & {:5.3f} & \\\\"
tabfoot="""\
    \end{{tabular}}
    \caption = "\caption{{Lines in the {reg} region, that extends from {start:.2f} to {finish:.2f} nm.}}\label{{tab:{reg}}}"
\end{{table}}"""

regions =["5053","5215","5654","6405","6449"]

for region in regions:
    reg = __import__(region)
    lfil = open("bin/{}_qu1.lin".format(region),"rb")
    lines = pic.load(lfil); lfil.close()
    print(tabhead)
    for line in lines:
        if line.El < 0:
            el = gf = "  -   "
        else:
            el = "{:5.3f} ".format(line.El)
            gf = "{:6.3f}".format(line.gf)
        print(tabline.format(line.name[:-7],line.cent,el,gf,line.dept))
    print(tabfoot.format(reg=region,start=reg.qu1.lmbd[-1],finish=reg.qu1.lmbd[0]))

