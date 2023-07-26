# import os
# os.environ["CUDA_VISIBLE_DEVICES"] = "0"

from __future__ import print_function
import argparse
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import torch
import torch.nn.parallel
import torch.optim as optim
import torch.utils.data
sys.path.append(os.getcwd())
# change 'E:\@\@\@\Pointnet.pytorch-master' below to your own project directory
sys.path.append('D:\\JASONZ\\FYP\\codepkg\\Pointnet')
from pointnet.dataset import ShapeNetDataset
from pointnet.model import PointNetCls, feature_transform_regularizer
from path2json import data_split, epoch_shuffle
from sklearn.metrics import confusion_matrix
from sklearn.metrics import ConfusionMatrixDisplay
import torch.nn.functional as F
from tqdm import tqdm

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--batchSize', type=int, default=3, help='input batch size')
    parser.add_argument(
        '--num_points', type=int, default=1000, help='input batch size')
    parser.add_argument(
        '--workers', type=int, help='number of data loading workers', default=6)
    parser.add_argument(
        '--nepoch', type=int, default=30, help='number of epochs to train for')
    parser.add_argument('--outf', type=str, default='cls', help='output folder')
    parser.add_argument('--model', type=str, default="D:\JASONZ\FYP\codepkg(2023.3.9)\pointnet_lstm\model70.pth", help='model path')
    parser.add_argument('--dataset', type=str, required=True, help="dataset path")
    parser.add_argument('--dataset_type', type=str, default='shapenet', help="dataset type shapenet|modelnet40")
    parser.add_argument('--feature_transform', action='store_true', help="use feature transform")
    blue = lambda x: '\033[94m' + x + '\033[0m'
    opt = parser.parse_args()
    print(opt)
    try:
        os.makedirs(opt.outf)
    except OSError:
        pass

    classifier = PointNetCls(k=4, feature_transform=opt.feature_transform)

    if opt.model != '':
        classifier.load_state_dict(torch.load(opt.model))

    if opt.dataset_type == 'shapenet':
        dataset = ShapeNetDataset(
            root=opt.dataset,
            classification=True,
            npoints=opt.num_points)
        test_dataset = ShapeNetDataset(
            root=opt.dataset,
            classification=True,
            split='test',
            npoints=opt.num_points,
            data_augmentation=False)
    else:
        exit('wrong dataset type')

    dataloader = torch.utils.data.DataLoader(
        dataset,
        batch_size=opt.batchSize,
        shuffle=True,
        drop_last=True,
        num_workers=int(opt.workers))

    testdataloader = torch.utils.data.DataLoader(
        test_dataset,
        drop_last=True,
        batch_size=opt.batchSize,
        shuffle=True,
        num_workers=int(opt.workers))
    
    optimizer = optim.Adam(classifier.parameters(), lr=0.0001, betas=(0.9, 0.999))
    scheduler = optim.lr_scheduler.StepLR(optimizer, step_size=10, gamma=0.5)

    classifier.cuda()
    total_correct = 0
    total_testset = 0
    target_list = []
    pred_list = []
    for i, data in tqdm(enumerate(testdataloader, 0)):
        points, target = data
        target = target[:, 0]
        target_list = np.append(target_list, target.cpu().numpy())
        points = points.transpose(2, 1)
        points, target = points.cuda(), target.cuda()
        classifier = classifier.eval()
        pred, _, _ = classifier(points)
        pred_choice = pred.data.max(1)[1]
        pred_list = np.append(pred_list, pred_choice.cpu().numpy())
        correct = pred_choice.eq(target.data).cpu().sum()
        total_correct += correct.item()
        total_testset += points.size()[0]
    cm = confusion_matrix(target_list, pred_list)
    print(cm)
    print("final accuracy {}".format(total_correct / float(total_testset)))
    lm = ['Jumping jack', 'Kick', 'Punch', 'Walk']
    disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=lm)
    disp.plot()
    plt.show()