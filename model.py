import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import Dataset
from torch.utils.data import DataLoader

class TransformerSetPartitionClassifier(nn.Module):
    def __init__(self, n, max_len, embedding_dim=32, num_heads=4, num_layers=2, hidden_dim=128):
        super().__init__()
        self.n = n
        self.embedding_dim = embedding_dim

        self.elem_embedding = nn.Embedding(n + 2, embedding_dim, padding_idx=0)  # +2 just to be safe
        self.block_embedding = nn.Embedding(max_len + 1, embedding_dim)
        self.r_embedding = nn.Linear(1, embedding_dim)

        encoder_layer = nn.TransformerEncoderLayer(d_model=embedding_dim * 3, nhead=num_heads, batch_first=True)
        self.s_transformer = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)

        j_encoder_layer = nn.TransformerEncoderLayer(d_model=embedding_dim, nhead=num_heads, batch_first=True)
        self.j_transformer = nn.TransformerEncoder(j_encoder_layer, num_layers=num_layers)

        self.output = nn.Linear(embedding_dim * 4, 3)

    def forward(self, S, r, J):
        device = next(self.parameters()).device  # Auto-move tensors to model device (cpu or cuda)

        s_tokens = []
        for block_id, (block, r_val) in enumerate(zip(S, r)):
            r_tensor = torch.tensor([r_val], dtype=torch.float, device=device).unsqueeze(0)
            r_emb = self.r_embedding(r_tensor).squeeze(0)
            b_emb = self.block_embedding(torch.tensor(block_id, device=device))

            for elem in block:
                e_emb = self.elem_embedding(torch.tensor(elem, device=device))
                token = torch.cat([e_emb, b_emb, r_emb], dim=0)
                s_tokens.append(token)

        s_tokens_tensor = torch.stack(s_tokens).unsqueeze(0)  # (1, len, dim)
        S_encoded = self.s_transformer(s_tokens_tensor)
        S_repr = S_encoded.mean(dim=1)  # (1, dim)

        j_tensor = torch.tensor(J, dtype=torch.long, device=device).unsqueeze(0)
        j_emb = self.elem_embedding(j_tensor)  # (1, len(J), dim)
        J_encoded = self.j_transformer(j_emb)
        J_repr = J_encoded.mean(dim=1)  # (1, dim)

        combined = torch.cat([S_repr, J_repr], dim=-1)  # (1, 2*dim)
        logits = self.output(combined.squeeze(0))  # (3,)
        return logits

    def forward_from_tensor(self, input_tensor, s_len, r_len, j_len):
        """
        Forward pass using a flat input tensor.
        """
        device = next(self.parameters()).device
        idx = 0

        # Recover S
        S_elems = input_tensor[idx:idx+s_len] * self.n
        idx += s_len

        # Recover r
        r_elems = input_tensor[idx:idx+r_len] * self.n
        idx += r_len

        # Recover J
        J_elems = input_tensor[idx:idx+j_len] * self.n

        # Round to integers
        S_list = [int(torch.round(x).item()) for x in S_elems]
        r_list = [int(torch.round(x).item()) for x in r_elems]
        J_list = [int(torch.round(x).item()) for x in J_elems]
        # Reconstruct the structure of S
        # Assume we know block sizes or use dummy blocks for now
        # (You can improve this later!)

        # Now use your normal forward(S, r, J)
        return self.forward(S=[S_list], r=r_list, J=J_list)

class SetPartitionDataset(Dataset):
    def __init__(self, data):
        self.data = data

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        (S, r, J), label = self.data[idx]
        return S, r, J, label