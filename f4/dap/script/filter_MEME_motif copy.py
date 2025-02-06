###run in python3

#goal: get the meme motif in the formate of Minimal MEME Motif Format (https://meme-suite.org/meme/doc/meme-format.html) for the FIMO input

#sys.argv[1]: the meme motif files that are the putput of the "MEME" running such as "/Users/ming/Dropbox/Kmer_manuscript/2021_NP/MEME_pipeline/MEME_results/G5_At_Discriminative/meme.txt"

#sys.argv[2]: the threshold for E-value (-evt) such as "5e-02"

#sys.argv[3]: meme_resource such "DIS" for Discovery Mode	Discriminative, "DE" for Discovery Mode	Differential Enrichment

##note_:
##this version for the fimo.tsv that is generated directly based on MEME meme.txt output


import sys

#input threshold of being a significant MEME motif
threshold_val=float(sys.argv[2])

#input threshold of being a significant MEME motif
meme_resource=(sys.argv[3])


in1=open(sys.argv[1], 'r')
out=open(sys.argv[1]+"_"+sys.argv[2]+"_"+sys.argv[3], 'w')

n=0
line_p=0
Print_back=False
ge_out=0

all_meme={}
all_meme_evalue={}
Print_motif=False

n_mo=0
motif_len=0
for l2 in in1:
    l1=l2.strip()
    #print ("/n"%(l1))
    n+=1
    if l1.startswith("MEME version"):
        print (l1)
        out.write("%s\n"%(l1))
    elif l1.startswith("ALPHABET") or l1.startswith("strands"):
        out.write("\n%s\n"%(l1))

    elif l1.startswith("Background letter frequencies"):
        out.write("\n%s\n"%(l1))
        Print_back=True
        line_p=n
    elif n == line_p+1 and Print_back:
        out.write("\n%s"%(l1))
        line_p=0
        Print_back=False

    elif l1.startswith("Motif") and "position-specific" in l1 and "probability" in l1:
        motif_se=""
        mot_in=l1.split()
        print (mot_in)

        motif_se=mot_in[2]+"_"+mot_in[1]+"_"+meme_resource
        all_meme[motif_se]=mot_in[1]
        n_mo=n
        motif_len=len(mot_in[1])
        motif_se=mot_in[2]+"_"+mot_in[1]+"_"+meme_resource


    elif n_mo < n < (n_mo+2+motif_len+1):
        if "letter-probability matrix" in l1 and float(l1.split()[-1]) < threshold_val:
            Print_motif=True
            out.write("\n\nMOTIF %s %s\n%s"%(motif_se, motif_se, l1))

            all_meme_evalue[motif_se]=(l1.split()[-1])

        elif Print_motif and "letter-probability matrix" not in l1:
            out.write("\n%s"%(l1))


    elif n > (n_mo+2+motif_len+1):
        Print_motif=False
out.write("\n")


#outfile_meme_list
out2=open(sys.argv[1]+"_"+sys.argv[2]+"_"+sys.argv[3]+"_motif_list", 'w')
#out.write("#python %s\n"%" ".join(sys.argv))
out2.write("meme_id\tmeme_Seq\tE_values\n")

for i_m, i_se in all_meme.items():
    if i_m in all_meme_evalue.keys():
    	out2.write("%s\t%s\t%s\n"%(i_m, i_se, all_meme_evalue[i_m]))



in1.close()
out.close()
out2.close()
