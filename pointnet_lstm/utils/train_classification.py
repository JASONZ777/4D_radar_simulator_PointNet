from __future__ import print_function
import argparse
import sys
import os
os.environ['CUDA_VISIBLE_DEVICES'] = '0'
import random
import torch
import torch.nn.parallel
import torch.optim as optim
import torch.utils.data
from torch.utils.tensorboard import SummaryWriter
sys.path.append(os.getcwd())
# change 'E:\@\@\@\Pointnet.pytorch-master' below to your own project directory
sys.path.append('D:\\JASONZ\\FYP\\codepkg(2023.3.9)\\pointnet_lstm')
from pointnet.dataset import ShapeNetDataset
from pointnet.model import PointNetCls, feature_transform_regularizer
from path2json import data_split, epoch_shuffle
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
        '--nepoch', type=int, default=200, help='number of epochs to train for')
    parser.add_argument('--outf', type=str, default='cls', help='output folder')
    parser.add_argument('--model', type=str, default='', help='model path')
    parser.add_argument('--dataset', type=str, required=True, help="dataset path")
    parser.add_argument('--dataset_type', type=str, default='shapenet', help="dataset type shapenet|modelnet40")
    parser.add_argument('--feature_transform', action='store_true', help="use feature transform")
    blue = lambda x: '\033[94m' + x + '\033[0m'
    opt = parser.parse_args()

    print(opt)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    print(device)
    opt.manualSeed = random.randint(1, 1000)  # fix seed
    print("Random Seed: ", opt.manualSeed)
    random.seed(opt.manualSeed)
    torch.manual_seed(opt.manualSeed)

    if opt.dataset_type == 'shapenet':
        dataset = ShapeNetDataset(
            root=opt.dataset,
            split='train',
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
    print(len(dataset), len(test_dataset))
    num_classes = len(dataset.classes)
    print('classes', num_classes)

    log = SummaryWriter()
    try:
        os.makedirs(opt.outf)
    except OSError:
        pass

    classifier = PointNetCls(k=4, feature_transform=opt.feature_transform)
    if opt.model != '':
        classifier.load_state_dict(torch.load(opt.model))

    classifier.cuda()
    optimizer = optim.Adam(classifier.parameters(), lr=0.0001, betas=(0.9, 0.999))
    scheduler = optim.lr_scheduler.StepLR(optimizer, step_size=10, gamma=0.5)
    num_batch = len(dataset) / opt.batchSize
    for epoch in range(opt.nepoch):
        #  scheduler.step()
        train_loss = 0
        train_right = 0
        test_loss = 0
        test_right = 0
        batch_test = 0
        for i, data in enumerate(dataloader, 0):
            points, target = data
            target = target[:, 0]
            points = points.transpose(2, 1)
            # points = points.to(device)
            # target = target.to(device)
            points, target = points.cuda(), target.cuda()
            optimizer.zero_grad()
            classifier = classifier.train()
            pred, trans, trans_feat = classifier(points)
            loss = F.nll_loss(pred, target)
            if opt.feature_transform:
                loss += feature_transform_regularizer(trans_feat) * 0.001
            loss.backward()
            optimizer.step()
            pred_choice = pred.data.max(1)[1]
            correct = pred_choice.eq(target.data).cpu().sum()
            print('[%d: %d/%d] train loss: %f accuracy: %f' % (
            epoch, i, num_batch, loss.item(), correct.item() / float(opt.batchSize)))
            train_loss += loss.item()   # calculate loss of each epoch
            train_right += correct.item()
            if i % 10 == 0:
                batch_test += 1
                j, data = next(enumerate(testdataloader, 0))
                points, target = data
                target = target[:, 0]
                points = points.transpose(2, 1)
                points, target = points.cuda(), target.cuda()
                classifier = classifier.eval()
                pred, _, _ = classifier(points)
                loss = F.nll_loss(pred, target)
                pred_choice = pred.data.max(1)[1]
                correct = pred_choice.eq(target.data).cpu().sum()
                print('[%d: %d/%d] %s loss: %f accuracy: %f' % (
                epoch, i, num_batch, blue('test'), loss.item(), correct.item() / float(opt.batchSize)))  
                test_loss += loss.item()  # calculate loss of each epoch
                test_right += correct.item()

        train_acc = train_right / len(dataset)
        test_acc = test_right / (batch_test*opt.batchSize)
        print(train_loss, train_right, test_loss, test_right)  
        info1 = {'train_loss': train_loss/num_batch, 'test_loss': test_loss/batch_test}
        log.add_scalars('Loss', info1, epoch)

        info2 = {'train_acc': train_acc, 'test_acc': test_acc}
        log.add_scalars('Accuracy', info2, epoch)   

        scheduler.step()
        if epoch % 10 == 0:
            torch.save(classifier.state_dict(), '%s\cls_model_%d.pth' % (opt.outf, epoch))

    total_correct = 0
    total_testset = 0
    for i, data in tqdm(enumerate(testdataloader, 0)):
        points, target = data
        target = target[:, 0]
        points = points.transpose(2, 1)
        points, target = points.cuda(), target.cuda()
        classifier = classifier.eval()
        pred, _, _ = classifier(points)
        pred_choice = pred.data.max(1)[1]
        correct = pred_choice.eq(target.data).cpu().sum()
        total_correct += correct.item()
        total_testset += points.size()[0]

    print("final accuracy {}".format(total_correct / float(total_testset)))