from captum.attr import IntegratedGradients
import training

def encode_input(S, r, J, max_n=10):
    """
    Encode (S, r, J) into a 1D tensor.
    max_n: maximum n we allow for padding
    """

    # Flatten S (list of lists) into a list of elements
    S_flat = [elem for block in S for elem in block]

    # Normalize everything to be between 0 and 1 (optional but helps)
    S_norm = [e / max_n for e in S_flat]
    r_norm = [float(val) / max_n for val in r]
    J_norm = [e / max_n for e in J]

    input_vec = S_norm + r_norm + J_norm
    input_tensor = torch.tensor(input_vec, dtype=torch.float32, requires_grad=True)

    return input_tensor

# Example input
S = [[1,3], [2,4,5]]
r = (1,1)
J = [1,3,5]

input_tensor = encode_input(S, r, J, max_n=5).to(device)
s_len = sum(len(block) for block in S)
r_len = len(r)
j_len = len(J)

ig = IntegratedGradients(model)

attr, delta = ig.attribute(
    inputs=input_tensor.unsqueeze(0),
    target=target_class_index,  # 0, 1, or 2
    additional_forward_args=(s_len, r_len, j_len),
    return_convergence_delta=True
)

s_attr = attr[0, :s_len]
r_attr = attr[0, s_len:s_len+r_len]
j_attr = attr[0, s_len+r_len:]

print("S attributions:", s_attr)
print("r attributions:", r_attr)
print("J attributions:", j_attr)