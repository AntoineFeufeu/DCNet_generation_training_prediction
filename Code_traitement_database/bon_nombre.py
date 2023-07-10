import os
# This code aims to putright numbers on a data base, so that all of the number samples in it starts with 0001 and follows each other.
stockage_data_name = "./ML_DATA/"
iend_database = 10000

for i in range (1,iend_database) :
    
    if not os.path.isdir(stockage_data_name+str(i).zfill(4)) :
        print(i)
        for k in range(i+1,iend_database) :
            if os.path.isdir(stockage_data_name+str(k).zfill(4)) :
                os.rename(stockage_data_name+str(k).zfill(4),stockage_data_name+str(i).zfill(4))
                break


# for i in range (2000,2718) :  
#     if os.path.isdir(stockage_data_name+str(i).zfill(4)) :
#         os.rename(stockage_data_name+str(i).zfill(4),stockage_data_name+str(2000+i).zfill(4))

