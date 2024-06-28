import torch
import torch.nn as nn
import torch.optim as optim
import math
from Bio.PDB import PDBParser, Select
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue

class RNAOptimizationModel(nn.Module):
    def __init__(self, n_units, pdb_file_path, spacing=500.0, device='cpu', alpha=1.0, lambda1=0.5, lambda2=0.5, delta=1.0, sigma=1.0, compression=False):
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
        self.euler_angles = nn.Parameter(torch.randn(n_units, 3, device=self.device))

        # Load coordinates from a PDB file
        self.parser = PDBParser(QUIET=True)
        self.structure = self.parser.get_structure('RNA', pdb_file_path)
        self.coordinates = []  # List to store coordinates
        self.extract_coordinates()
        if compression:
            self.compress_coordinates()

    def extract_coordinates(self):
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    if 'P' in residue and 'C1\'' in residue:
                        p_coord = residue['P'].get_vector().tolist()
                        c1_coord = residue['C1\''].get_vector().tolist()
                        self.coordinates.extend([p_coord, c1_coord])
        self.coordinates = torch.tensor(self.coordinates, dtype=torch.float, device=self.device)

    def compress_coordinates(self):
        # Placeholder: Implement actual compression logic if needed
        pass

    def forward(self):
        rotations = self.compute_rotation_matrices(self.euler_angles)
        loss = self.compute_loss(self.positions, rotations)
        return loss

    def compute_rotation_matrices(self, euler_angles):
        phi, theta, psi = euler_angles[:, 0], euler_angles[:, 1], euler_angles[:, 2]
        cos_phi, sin_phi = torch.cos(phi), torch.sin(phi)
        cos_theta, sin_theta = torch.cos(theta), torch.sin(theta)
        cos_psi, sin_psi = torch.cos(psi), torch.sin(psi)
        Rz = torch.stack([torch.stack([cos_psi, -sin_psi, torch.zeros_like(psi)], dim=-1),
                          torch.stack([sin_psi, cos_psi, torch.zeros_like(psi)], dim=-1),
                          torch.stack([torch.zeros_like(psi), torch.zeros_like(psi), torch.ones_like(psi)], dim=-1)], dim=-2)
        Ry = torch.stack([torch.stack([cos_theta, torch.zeros_like(theta), sin_theta], dim=-1),
                          torch.stack([torch.zeros_like(theta), torch.ones_like(theta), torch.zeros_like(theta)], dim=-1),
                          torch.stack([-sin_theta, torch.zeros_like(theta), cos_theta], dim=-1)], dim=-2)
        Rx = torch.stack([torch.stack([torch.ones_like(phi), torch.zeros_like(phi), torch.zeros_like(phi)], dim=-1),
                          torch.stack([torch.zeros_like(phi), cos_phi, -sin_phi], dim=-1),
                          torch.stack([torch.zeros_like(phi), sin_phi, cos_phi], dim=-1)], dim=-2)
        R = torch.matmul(Rz, Ry)
        R = torch.matmul(R, Rx)
        return R

    def compute_loss(self, positions, rotations):
        transformed = torch.matmul(rotations, self.coordinates.unsqueeze(-1)).squeeze(-1) + positions.unsqueeze(1)
        # L1 term
        L1 = torch.mean(torch.sum(transformed**2, dim=[1, 2]))
        # Global center of mass term
        global_com = torch.sum(self.positions, dim=0)
        L_center_of_mass = self.alpha * torch.sum(global_com**2)
        # Interaction terms (not fully implemented here, placeholder for now)
        L_interaction = torch.tensor(0.0, device=self.device)  # Placeholder for interaction term calculations
        return L1 + L_center_of_mass + L_interaction
    def save_pdb(self, output_path):
        # Set up PDB structure
        io = PDBIO()
        s = self.structure[0]  # Use the first model in the loaded structure as a template
        new_structure = s.copy()
        new_structure.detach_parent()
        
        # Generate new chains based on the optimized positions and rotations
        chain_id = 'A'
        for i, (pos, rot) in enumerate(zip(self.positions, self.compute_rotation_matrices(self.euler_angles))):
            chain = Chain(chain_id)
            new_structure.add(chain)
            for j, residue in enumerate(s.get_residues()):
                new_residue = Residue(residue.id, residue.resname, residue.segid)
                chain.add(new_residue)
                for atom in residue:
                    new_atom = atom.copy()
                    new_atom.transform(rot, pos)
                    new_residue.add(new_atom)
            
            chain_id = chr(ord(chain_id) + 1)  # Increment chain ID
            
        # Save the new structure to a PDB file
        io.set_structure(new_structure)
        io.save(output_path)

        print(f"Saved optimized structure to {output_path}")

# Device setup
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
# Model instantiation
input_monomer_pdb = "/large/otgk/casp/casp16/2_R1205/pdb/S_000001.pdb"
n = 3
model = RNAOptimizationModel(n_units=n, pdb_file_path=input_monomer_pdb, device=device)
# Optimizer setup
optimizer = optim.Adam(model.parameters(), lr=0.01)
# Optimization loop
for epoch in range(1000):
    print(f"Epoch {epoch}")
    optimizer.zero_grad()
    loss = model()
    loss.backward()
    optimizer.step()
    if epoch % 100 == 0:
        print(f"Epoch {epoch}: Loss = {loss.item()}")
# Save optimized structure to a PDB file
output_pdb = f"/large/otgk/casp/casp16/2_R1205/pdb/multimer/S_000001_optimized.{n}mer.pdb"
parameter_file = f"/large/otgk/casp/casp16/2_R1205/pdb/multimer/S_000001_optimized.{n}mer.params"
model.save_pdb(output_pdb)
# Save parameters to a file
torch.save(model.state_dict(), parameter_file)
