from Bio import SeqIO

flag_kos = ["K02387","K02388","K02389","K02390","K02391","K02392","K02393","K02394","K02396","K02397","K02400","K02401","K02406","K02407","K02408","K02409","K02410","K02411","K02412","K02414","K02416","K02417","K02418","K02419","K02420","K02421","K02556","K02557","K13820","K21217","K21218","K02386","K02399","K02413","K02422","K02423"]

with open("wsip/results/dram/metawrap-drep-bins/genbank.files.txt") as file:
    gbks = file.readlines()
    gbks = [line.rstrip() for line in gbks]

for gbk in gbks:
    for record in SeqIO.parse("wsip/results/dram/metawrap-drep-bins/genbank/" + gbk,"genbank"):
        for feature in record.features:
            for ko in flag_kos:
                myko = "kegg:" + ko
                if feature.type == "CDS" and myko in feature.qualifiers.get("db_xref",[]):
                    printline = gbk + "," + ko + "," + feature.qualifiers.get("gene")[0]
                    print(printline)
