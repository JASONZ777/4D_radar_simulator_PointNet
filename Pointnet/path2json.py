import os
import json
import math
import random
# 打开文件


def data_split(windowsize, slidingsize):
    lis = [[], [], []]
    newlist = []
    traingroup = []
    testgroup = []
    split_train = []
    split_test = []
    classlist = ['jumpingJack3Reps/', 'kickLSide2Reps/', 'walkLeft3Steps/']
    for i in range(0, len(classlist)):    # 将不同活动文件读取并放到同一个file里面     
        COOKED_FOLDER = 'acty_add_velocity/' + classlist[i]
        dirs = os.listdir(COOKED_FOLDER)
        dirs.sort(key=lambda x: (int(x.split('.')[0].split('_')[1]), int(x.split('.')[0].split('_')[2])))  # 排序！！！python读取文件夹文件可能会乱序
        # 输出所有文件和文件夹
        for file in dirs:
            filepath = COOKED_FOLDER + file
            lis[i].append(filepath)

        n_newlist = range(0, math.floor((len(lis[i])-windowsize)/slidingsize)+1, 1)  # 设定时间序列窗大小，这里设置滑动overlap

        for n in n_newlist:
            for t in range(slidingsize*n, slidingsize*n+windowsize, 1):
                newlist.append(lis[i][t])
    group = [newlist[i: i+windowsize] for i in range(0, len(newlist), windowsize)]  # 设置 group，为了能够shuffle
    manualSeed = random.randint(1, 1000)  # fix seed
    random.seed(manualSeed)
    random.shuffle(group)
    for n in range(0, len(group), 1):   # 划分训练集和测试集
        if n % 10 != 0:
            traingroup.extend(group[n])  # 追加，用于写入
            split_train.append(group[n])  # 扩展， 用于shuffle
        else:
            testgroup.extend(group[n])
            split_test.append(group[n])
     
    with open('acty_add_velocity/train_file.json', 'w') as f:   # 存入训练集 9：1
        json.dump(traingroup, f)
    with open('acty_add_velocity/test_file.json', 'w') as f:   # 存入测试集
        json.dump(testgroup, f)
    return split_train, split_test


def epoch_shuffle(split_train, split_test):
    manualSeed = random.randint(1, 1000)  # fix seed
    print("Random Seed: ", manualSeed)
    random.seed(manualSeed)
    random.shuffle(split_train)
    random.shuffle(split_test)
    split_test = [i for item in split_test for i in item]
    split_train = [i for item in split_train for i in item]

    with open('acty_add_velocity/train_file.json', 'w') as f:   # 存入训练集 9：1
        json.dump(split_train, f)
    with open('acty_add_velocity/test_file.json', 'w') as f:   # 存入测试集
        json.dump(split_test, f)


