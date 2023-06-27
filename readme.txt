#0. Install Python libraries needed

###The python environment is exact same as GraphDTA(https://github.com/thinng/GraphDTA)

conda activate geometric

1. prepare the data
####the original data are in positive_n and negative_n

python read_smi_protein_nnn_pos.py

2. run the prediction
##unzip full_model_out2000.rar file and put the full_model_out2000.model in the main folder.
python training_nn3_load_name.py


####nohup.out_MCC is the log file, and output_f.txt is the collected information by runing "python  read_out.py"
This work3_VS_n_general_O14626_gcn_BC_cutoff0.6_upload.zip is a example of screening chemdiv dataase for O14626, but we only provided few compound smiles, beause whole chemdiv is too large to deposite in github. User can prepare their all dataset using this code for large scale virtual screening.
