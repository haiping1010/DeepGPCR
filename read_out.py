
fr=open('nohup.out_MCC','r')

arr=fr.readlines()

print ('AUC:'+'\t'+'TPR:'+'\t'+'precision'+'\t'+'accuracy'+'\t'+'MCC'+'\t'+'pos_num'+'\t'+'neg_num')
for line in arr:
    if  line.startswith('AUC:'):
        tem_arr=line.split(':')
        out_line=tem_arr[1].strip()
    if  line.startswith('TPR:'):
        tem_arr=line.split(':')
        out_line=out_line+'\t'+tem_arr[1].strip()
    if  line.startswith('precision:'):
        tem_arr=line.split()
        out_line=out_line+'\t'+tem_arr[1].strip()
    if  line.startswith('accuracy:'):
        tem_arr=line.split()
        out_line=out_line+'\t'+tem_arr[1].strip()
    if  line.startswith('MCC:'):
        tem_arr=line.split()
        out_line=out_line+'\t'+tem_arr[1].strip()


    if 'pos_num:'  in line:
        arr_tem=line.split(':')
        pos=arr_tem[1].strip()
        out_line=out_line+'\t'+pos
    if 'neg_num:'  in line:
        arr_tem=line.split(':')
        neg=arr_tem[1].strip()
        out_line=out_line+'\t'+neg

        print (out_line.strip())




