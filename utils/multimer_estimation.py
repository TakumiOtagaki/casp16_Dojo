import matplotlib.pyplot as plt
import torch
import torch.nn as nn
import torch.optim as optim
import math
from Bio.PDB import PDBParser, Select
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import Chain
from Bio.PDB import Residue
from Bio.PDB import Atom
from Bio.PDB import Model
from Bio.PDB import Structure
import numpy as np
import sys
torch.autograd.set_detect_anomaly(True)


class RNAOptimizationModel(nn.Module):
    def __init__(self, n_units, pdb_file_path, spacing=300.0, device='cpu', alpha=1.0, lambda1=10, lambda2=10, delta=1.0, sigma=1.0, compression=True):
        super(RNAOptimizationModel, self).__init__()
        self.n_units = n_units
        self.device = device
        self.alpha = alpha
        self.lambda1 = lambda1
        self.lambda2 = lambda2
        self.delta = delta
        self.sigma = sigma

        # Initialize positions in a circular layout
        positions = torch.zeros(n_units, 3, device=device)
        angle = 2 * math.pi / n_units
        for i in range(n_units):
            positions[i, 0] = spacing * math.cos(i * angle)
            positions[i, 1] = spacing * math.sin(i * angle)
            positions[i, 2] = torch.randn(1, device=device) * 0.1
        self.positions = nn.Parameter(positions)

        # Initialize Euler angles for rotation matrices
        self.euler_angles = nn.Parameter(
            torch.randn(n_units, 3, device=self.device))

        # Load coordinates from a PDB file
        self.parser = PDBParser(QUIET=True)
        self.structure = self.parser.get_structure('RNA', pdb_file_path)
        self.coordinates = []  # List to store coordinates
        self.extract_coordinates()
        print(f"Extracted {len(self.coordinates)} coordinates from PDB file")
        if compression:
            self.compress_coordinates()
        print(f"Compressed coordinates to {len(self.coordinates)} atoms")
        print("box start")
        self.calculate_boxes()

    def extract_coordinates(self):
        # valid_residues =
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    # if 'P' in residue and 'C1\'' in residue:
                    #     p_coord = residue['P'].get_vector()
                    #     c1_coord = residue['C1\''].get_vector()
                    #     p_coord = torch.tensor(
                    #         [p_coord[0], p_coord[1], p_coord[2]], device=self.device)
                    #     c1_coord = torch.tensor(
                    #         [c1_coord[0], c1_coord[1], c1_coord[2]], device=self.device)
                    #     self.coordinates.extend([p_coord, c1_coord])
                    # C1 だけにする
                    if 'C1\'' in residue:
                        c1_coord = residue['C1\''].get_vector()
                        c1_coord = torch.tensor(
                            [c1_coord[0], c1_coord[1], c1_coord[2]], device=self.device)
                        self.coordinates.append(c1_coord)
        self.coordinates = torch.stack(self.coordinates).view(-1, 3).float()

    def compress_coordinates(self, threshold=0.70):
        initial_count = len(self.coordinates)
        compressed_coords = []
        i = 0
        while i < len(self.coordinates) - 2:
            u = self.coordinates[i + 1] - self.coordinates[i]
            v = self.coordinates[i + 2] - self.coordinates[i + 1]
            u_norm = u / torch.norm(u)
            v_norm = v / torch.norm(v)
            cos_similarity = torch.dot(u_norm, v_norm)

            if cos_similarity > threshold:
                i += 1  # Skip the middle point
            else:
                compressed_coords.append(self.coordinates[i])
            i += 1

        # Add the last two coordinates if not already added
        compressed_coords.extend(self.coordinates[-2:])
        self.coordinates = torch.stack(compressed_coords)
        print(
            f"Compressed from {initial_count} to {len(self.coordinates)} atoms")

    def calculate_boxes(self):
        print("hi")
        min_coords = torch.min(self.coordinates, dim=0).values
        max_coords = torch.max(self.coordinates, dim=0).values
        self.box = torch.stack([min_coords, max_coords], dim=0)  # [2, 3]
        print("end hi")

    def forward(self):
        rotations = self.compute_rotation_matrices(self.euler_angles)
        loss = self.compute_loss(self.positions, rotations)
        return loss

    def compute_rotation_matrices(self, euler_angles):
        phi, theta, psi = euler_angles[:,
                                       0], euler_angles[:, 1], euler_angles[:, 2]
        cos_phi, sin_phi = torch.cos(phi), torch.sin(phi)
        cos_theta, sin_theta = torch.cos(theta), torch.sin(theta)
        cos_psi, sin_psi = torch.cos(psi), torch.sin(psi)
        Rz = torch.stack([torch.stack([cos_psi, -sin_psi, torch.zeros_like(psi)], dim=-1),
                          torch.stack(
            [sin_psi, cos_psi, torch.zeros_like(psi)], dim=-1),
            torch.stack([torch.zeros_like(psi), torch.zeros_like(psi), torch.ones_like(psi)], dim=-1)], dim=-2)
        Ry = torch.stack([torch.stack([cos_theta, torch.zeros_like(theta), sin_theta], dim=-1),
                          torch.stack([torch.zeros_like(theta), torch.ones_like(
                              theta), torch.zeros_like(theta)], dim=-1),
                          torch.stack([-sin_theta, torch.zeros_like(theta), cos_theta], dim=-1)], dim=-2)
        Rx = torch.stack([torch.stack([torch.ones_like(phi), torch.zeros_like(phi), torch.zeros_like(phi)], dim=-1),
                          torch.stack(
            [torch.zeros_like(phi), cos_phi, -sin_phi], dim=-1),
            torch.stack([torch.zeros_like(phi), sin_phi, cos_phi], dim=-1)], dim=-2)
        R = torch.matmul(Rz, Ry)
        R = torch.matmul(R, Rx)
        return R

    def compute_loss(self, positions, rotations):
        m = self.coordinates.size(0)  # Number of coordinates
        n = self.positions.size(0)  # Number of units
        coordinates_expanded = self.coordinates.unsqueeze(0).repeat(
            rotations.size(0), 1, 1)  # [n_units, n_coordinates, 3]

        # 回転を適用
        transformed = torch.matmul(rotations, coordinates_expanded.transpose(
            1, 2)).transpose(1, 2)  # [n_units, n_coordinates, 3]

        # 位置ベクトルを追加するために形状を合わせる
        positions_expanded = positions.unsqueeze(1)  # [n_units, 1, 3]
        transformed = transformed + positions_expanded  # [n_units, n_coordinates, 3]
        # L1 term
        L1 = torch.mean(torch.sum(transformed**2, dim=[1, 2]))
        # Global center of mass term
        global_com = torch.sum(self.positions, dim=0)
        L_center_of_mass = self.alpha * torch.sum(global_com**2)
        # Interaction terms (not fully implemented here, placeholder for now)

        # Placeholder for interaction term calculations
        L_center_of_mass = self.alpha * torch.sum(global_com**2)

        # Interaction term
        print("L_interaction")
        L_interaction = torch.tensor(0.0, device=self.device)
        for i in range(n):
            for j in range(i + 1, n):
                for x1 in range(m):
                    for x2 in range(m):
                        dist_squared = torch.sum(
                            (transformed[i, x1] - transformed[j, x2])**2)
                        # Gaussian penalty for being too close
                        penalty = torch.exp(-(dist_squared -
                                              self.delta**2) / self.sigma**2)
                        # L_interaction += self.lambda1 * penalty
                        L_interaction = L_interaction + self.lambda1 * penalty

        # L_interaction / (n * (n - 1) * m**2)
        L_interaction = L_interaction / (n * (n - 1) * m**2)
        print("L_overlap")
        L_overlap = self.calculate_overlap_loss(
            positions, rotations) * self.lambda2 * 3

        total_loss = L1 + L_center_of_mass + L_interaction + L_overlap

        return total_loss

    def calculate_overlap_loss(self, positions, rotations):
        n = positions.size(0)
        L_overlap = torch.tensor(0.0, device=self.device)
        boxes = []

        for i in range(n):
            box = self.box.clone()
            # box[0] = torch.matmul(rotations[i], box[0]) + positions[i]
            # box[1] = torch.matmul(rotations[i], box[1]) + positions[i]
            tmp = torch.zeros(2, 3, device=self.device)
            tmp[0] = tmp[0] + torch.matmul(rotations[i], box[0]) + positions[i]
            tmp[1] = tmp[1] + torch.matmul(rotations[i], box[1]) + positions[i]
            box = box + tmp
            
            boxes.append(box)

        for i in range(n):
            for j in range(i + 1, n):
                overlap = self.box_overlap_volume(boxes[i], boxes[j])
                L_overlap = L_overlap + overlap

        return L_overlap

    def box_overlap_volume(self, box1, box2):
        min_corner1, max_corner1 = box1
        min_corner2, max_corner2 = box2

        overlap_min = torch.max(min_corner1, min_corner2)
        overlap_max = torch.min(max_corner1, max_corner2)
        overlap_dims = torch.clamp(overlap_max - overlap_min, min=0)
        overlap_volume = torch.prod(overlap_dims)

        return overlap_volume
# from Bio.PDB import PDBIO, Chain, Residue, Atom, Model, Structure

    def save_pdb(self, output_path):
        io = PDBIO()
        # 新しい空の構造を作成
        new_structure = Structure.Structure('optimized_structure')
        model = Model.Model(0)  # 新しいモデルを作成
        new_structure.add(model)

        # 利用可能なChain IDを生成
        available_chain_ids = (id for id in map(
            chr, range(65, 91)))  # 'A' to 'Z'

        # 回転行列と位置ベクトルを用いて各chainの座標を更新
        for i, (position, rotation) in enumerate(zip(self.positions.detach().cpu().numpy(), self.compute_rotation_matrices(self.euler_angles).detach().cpu().numpy())):
            chain_id = next(available_chain_ids, None)
            if chain_id is None:
                raise ValueError("Available chain IDs exhausted.")
            chain = Chain.Chain(chain_id)  # 新しいChainインスタンスを作成
            model.add(chain)

            # チェーンの各残基と原子に対して変換を適用
            for residue_id, residue in enumerate(self.structure.get_residues(), start=1):
                new_residue = Residue.Residue(
                    (' ', residue_id, ' '), residue.get_resname(), residue.segid)
                chain.add(new_residue)
                for atom in residue:
                    # 各原子に対して回転と平行移動を適用
                    transformed_coord = np.dot(rotation, atom.coord) + position
                    new_atom = Atom.Atom(atom.get_name(), transformed_coord, atom.get_bfactor(
                    ), atom.get_occupancy(), atom.get_altloc(), atom.get_fullname(), atom.get_serial_number())
                    new_residue.add(new_atom)

        # 新しい構造をPDBファイルに保存
        io.set_structure(new_structure)
        io.save(output_path)
        print(f"Saved optimized structure to {output_path}")


def optimize_and_plot(n, epochs, input_monomer_pdb, output_dir, device):
    model = RNAOptimizationModel(
        n_units=n, pdb_file_path=input_monomer_pdb, device=device)
    optimizer = optim.Adam(model.parameters(), lr=15)
    losses = []

    for epoch in range(epochs):
        # 初めは lr を大きくしておく
        if epoch == 100:
            for param_group in optimizer.param_groups:
                param_group['lr'] = 0.1
        elif epoch == 800:
            for param_group in optimizer.param_groups:
                param_group['lr'] = 0.01
        optimizer.zero_grad()
        loss = model()
        loss.backward()
        optimizer.step()
        losses.append(loss.item())
        print(f"N: {n}, Epoch {epoch + 1}/{epochs}, Loss: {loss.item()}")
        # # 定期的に pdb を保存
        # if (epoch + 1) % 100 == 0:
        if epoch % 100 == 0:
            output_pdb = f"{output_dir}/S_000001_optimized.{n}mer_epoch{epoch}.pdb"
            model.save_pdb(output_pdb)

    # Save the model and parameters
    output_pdb = f"{output_dir}/S_000001_optimized.{n}mer.pdb"
    parameter_file = f"{output_dir}/S_000001_optimized.{n}mer.params"
    model.save_pdb(output_pdb)
    torch.save(model.state_dict(), parameter_file)

        # Store losses for plotting
    return losses

# Parameters
n_values = range(2, 11)
epochs = 1000
# pjt_dir = "."
pjt_dir = "/work/gs58/s58007/casp16_Dojo" # wisteria
input_monomer_pdb = pjt_dir + "/R0254/pdb/S_000144.pdb"
output_dir = pjt_dir + "/R0254/pdb/multimer"

# Device setup
# device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
device = torch.device("cpu")

# Run optimization and plotting
import multiprocessing as mp
with mp.Pool(9) as pool:
    results = pool.starmap(optimize_and_plot, [(n, epochs, input_monomer_pdb, output_dir, device) for n in n_values])

# Plot the losses
plt.figure(figsize=(10, 6))
for n, losses in zip(n_values, results):
    plt.plot(losses, label=f"{n}-mer")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.yscale("log")
plt.legend()
plt.title("Losses during optimization")
plt.savefig(f"{output_dir}/losses.png")

