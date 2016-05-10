import pickle as pic

sectext = r"\subsection{{ The {} region}}"
figtext = """
\\begin{{figure}}[Th]
        \centering
        %trim option's parameter order: left bottom right top
        \includegraphics[trim = 6mm 2mm 2mm 2mm, clip,width=0.8\\textwidth]{{{reg}_{linname}_qu1.png}}
                \caption{{ {line} {wave} line }}
                \label{{ fig:{reg}qu1{linname} }}
\end{{figure}}"""

regnames = ["5053","5215","5654","6405","6449"]

for region in regnames:
    lfil = open("bin/{}_qu1.lin".format(region),"rb")
    lines = pic.load(lfil); lfil.close()
    print("")
    print(sectext.format(region))
    for line in lines:
        linname = "{}_{}".format(line.name.split()[0],int(line.cent*10))        
        name = line.name[:-7].strip()
        wave = int(line.cent*10)/10
        props = {"reg":region,"linname":linname,"line":name,"wave":wave}
        print(figtext.format(**props))
