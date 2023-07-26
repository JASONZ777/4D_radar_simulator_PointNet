import os
import numpy as np
import json
import pandas as pd
from sklearn.model_selection import KFold
# 打开文件

lis = []
train_files = []
test_files = []
classlist = ['jumpingJack3Reps/', 'kickLSide2Reps/', 'punchLSide2Reps/', 'walkLeft3Steps/']
for i in range(0, len(classlist)):    # 将不同活动文件读取并放到同一个file里面     
    COOKED_FOLDER = 'actyclass/' + classlist[i]
    dirs = os.listdir(COOKED_FOLDER)
    # 输出所有文件和文件夹
    for file in dirs:
        filepath = COOKED_FOLDER + file
        lis.append(filepath)
kf = KFold(n_splits=10, shuffle=True, random_state=42)  # 初始化KFold 10折
for k, (Trindex, Tsindex) in enumerate(kf.split(lis)):       
    train_files.append(np.array(lis)[Trindex].tolist())
    test_files.append(np.array(lis)[Tsindex].tolist())
print(len(train_files[0]))

with open('actyclass/train_file.json', 'w') as f:   # 存入训练集 9：1
    json.dump(train_files[0], f)
with open('actyclass/test_file.json', 'w') as f:   # 存入测试集
    json.dump(test_files[0], f)

