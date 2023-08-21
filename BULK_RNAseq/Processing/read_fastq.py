filename="../../../../Glioma_IGH/Bulk_Gliome/LGG85_GF/new_LGG85_GF_1.fq"

cts=0
with open(filename) as f:
    one_file = f.readlines()
    print(" Nb of lines : " + str(len(one_file)))
    for line in one_file:
        cts+=1
        if( cts%4==1):
            if not ( line.startswith('@')):
                print(cts)
                print(line)
        if(cts >  87425250 and cts <  87425270):
            print(line)

# head -n-2 LGG85_GF_1.fq > new_LGG85_GF_1.fq