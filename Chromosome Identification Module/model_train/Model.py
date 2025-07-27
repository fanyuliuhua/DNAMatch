import torch
import torch.nn.functional as F
import torch.nn as nn
from einops import rearrange, reduce
# from thop import profile


class Residual(nn.Module):
    def __init__(self, fn):
        super().__init__()
        self.fn = fn

    def forward(self, x, **kwargs):
        return self.fn(x, **kwargs) + x
class LinearNet(nn.Module):
    def __init__(self):
        super(LinearNet, self).__init__()
        # depth=5
        LNet=[]
        layer1=nn.Sequential(nn.Linear(768,3072),
                            nn.LeakyReLU(0.1),
                            nn.Dropout(0.5),
                            nn.Linear(3072,3072),
                             nn.LeakyReLU(0.1),
                             nn.Dropout(0.5),
                            nn.Linear(3072,768),
                            nn.LeakyReLU(),
                            )
        layer2=nn.Sequential(nn.Linear(768,3072),
                            nn.LeakyReLU(0.1),
                            nn.Dropout(0.4),
                            nn.Linear(3072,3072),
                             nn.LeakyReLU(0.1),
                             nn.Dropout(0.4),
                            nn.Linear(3072,768),
                            nn.LeakyReLU(),
                            )
        layer3=nn.Sequential(nn.Linear(768,3072),
                            nn.LeakyReLU(0.1),
                            nn.Dropout(0.3),
                            nn.Linear(3072,3072),
                             nn.LeakyReLU(0.1),
                             nn.Dropout(0.3),
                            nn.Linear(3072,768),
                            nn.LeakyReLU(),
                            )
        layer4=nn.Sequential(nn.Linear(768,3072),
                            nn.LeakyReLU(0.1),
                            nn.Dropout(0.2),
                            nn.Linear(3072,3072),
                             nn.LeakyReLU(0.1),
                             nn.Dropout(0.2),
                            nn.Linear(3072,768),
                            nn.LeakyReLU(),
                            )
        layer5=nn.Sequential(nn.Linear(768,3072),
                            nn.LeakyReLU(0.1),
                            nn.Dropout(0.1),
                            nn.Linear(3072,3072),
                             nn.LeakyReLU(0.1),
                             nn.Dropout(0.1),
                            nn.Linear(3072,768),
                            nn.LeakyReLU(),
                            )
        LNet.append(nn.Sequential(
                nn.BatchNorm1d(768),
                Residual(
                    layer1,
                         ),
                nn.Linear(768,768),
                nn.LeakyReLU(),
            nn.BatchNorm1d(768),
            Residual(
                layer2,
            ),
            nn.Linear(768, 768),
            nn.LeakyReLU(),
            nn.BatchNorm1d(768),
            Residual(
                layer3,
            ),
            nn.Linear(768, 768),
            nn.LeakyReLU(),
            nn.BatchNorm1d(768),
            Residual(
                layer4,
            ),
            nn.Linear(768, 768),
            nn.LeakyReLU(),
            nn.BatchNorm1d(768),
            Residual(
                layer5,
            ),
            nn.Linear(768, 768),
            nn.LeakyReLU(),
            ))
        self.LNet=nn.Sequential(*LNet)
        self.FinalLayer=nn.Sequential(
            nn.Linear(768,512),
            nn.Dropout(0.2),
            nn.LeakyReLU(),
            nn.Linear(512,256),
            nn.Dropout(0.1),
            nn.LeakyReLU(),
            nn.Linear(256,24)
        )
    def forward(self,x):
        x=self.LNet(x)
        x=self.FinalLayer(x)
        return x

# model=LinearNet()
# inputs=torch.randn((1000,768)).type(dtype=torch.float32)
# flops,params=profile(model,inputs=(inputs,))
# print(params)
# print(flops)