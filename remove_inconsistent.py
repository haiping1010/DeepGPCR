
import glob
import os
arr_pos=glob.glob("positive_n/??????/*.sdf")

arr_neg=glob.glob("negative_n/??????/*.sdf")

fw=open("removed_inconsistented_pos_neg.txt",'w')
os.system("mkdir  inconsisted")
for pos_name in arr_pos:
    for  neg_name  in arr_neg:
         if  pos_name.replace('positive_n/','') == neg_name.replace('negative_n/',''):
             fw.write(pos_name+"\t"+neg_name)
             fw.write("\n")
             print (pos_name,neg_name)
             os.system("rm -rf  "+pos_name.replace('.sdf','')+"*")
             os.system("rm -rf  "+neg_name.replace('.sdf','')+"*")




fw.close()
             
            

