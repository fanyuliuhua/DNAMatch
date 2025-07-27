import argparse
import torch



parser=argparse.ArgumentParser()


parser.add_argument('--d-model', default=768, type=int)
parser.add_argument('--gpu', default=0, type=int)
parser.add_argument('--lr', type=float, default=0.00001, help='Initial learning rate.')
parser.add_argument('--weight_decay', type=float, default=1e-7, help='Weight decay (L2 loss on parameters).')
args = parser.parse_args()
device = torch.device(f"cuda:{args.gpu}" if torch.cuda.is_available() else "cpu")
args.device = device