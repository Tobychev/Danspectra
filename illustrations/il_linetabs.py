import pickle as pic

tabhead="""
\\begin{table}[]
\centering
    \\begin{tabular}{lccccl}
Species & $\lambda$  & E${}_l$ & $\log gf$ & depth & Comment\\\\ \hline"""
tabline = "{}    & {:7.3f}    & {:5.3f}   & {:6.3f}    & {:5.3f} & \\\\"
tabfoot="""\
    \end{{tabular}}
    \caption{{Lines in the {reg} region.}}\label{{tab:{reg}}}
\end{{table}}"""

regions =["5053","5215","5654","6405","6449"]

for region in regions:
    lfil = open("bin/{}_qu1.lin".format(region),"rb")
    lines = pic.load(lfil); lfil.close()
    print(tabhead)
    for line in lines:
        print(tabline.format(line.name[:-7],line.cent,line.El,line.gf,line.dept))
    print(tabfoot.format(reg=region))

