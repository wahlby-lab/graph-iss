""" Signal candidate probability predictions.
    
    Parameters
    ----------
    sys.argv[1] : Path to signal probability prediction library
    sys.argv[2] : Input dataframe of signal candidate detection
        after merging
    sys.argv[3] : Output dataframe of signal probability predictions
    sys.argv[4] : (Optional) probability threshold for signal candid-
	date predictions
"""

from keras.models import model_from_json
import pandas as pd
import numpy as np
import sys


#loading data to classify and to be classified  
path=sys.argv[1]
model_name='/model_FT3'
try:
    m=pd.read_pickle(sys.argv[2])
except EOFError:
    sys.exit(0)

x_data=m['Intensities_window_5x5'] 
for i in range(len(x_data)):
    if x_data[i].size!=25 or np.isnan(x_data[i]).any(): #avoid window of intensities belonging to the borders
       x_data=x_data.drop([i])
       m=m.drop([i])
x_data=x_data.reset_index(drop=True)
m=m.reset_index(drop=True)

X_data=np.zeros([len(x_data),5,5])

for i in range(len(x_data)):
        X_data[i,:,:]=x_data[i] 
        
X_data=X_data.reshape(len(X_data),5,5,1)  

#getting probabilities from the DNN
json_file = open(path+model_name+'.json', 'r')
loaded_model_json = json_file.read()
json_file.close()
DNN_model = model_from_json(loaded_model_json)
DNN_model.load_weights(path+model_name+'_weights')
y_pred_DNN=DNN_model.predict_proba(X_data)

prob_df=pd.DataFrame(data={'prob_DNN':y_pred_DNN[:,1]})
final_df=pd.concat([m, prob_df],axis=1)

#saving relative signals probabilities
if len(sys.argv)>4:
    th = float(sys.argv[4])
    final_df = final_df[final_df.prob_DNN > th]
final_df.drop_duplicates(['cycle','ch','x','y','z'], inplace=True)
final_df = final_df.reset_index(drop=True)
final_df.to_pickle(sys.argv[3])
