#!/usr/bin/env python3
from random import sample
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt 
from tqdm import tqdm
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report
from sklearn.decomposition import PCA

from comet import mute_rdkit_warnings, smiles_to_mol, mol_to_barcode, bit_enrichment


def main() -> None:
    """
    Driver code for script.
    """
    mute_rdkit_warnings()

    # Hyperparameters
    test_size = 0.5
    random_state = 42 
    num_bits = 2084
    radius = 2
    alpha = 0.025
    mtc = "bonferroni"

    # Parase data from input file
    data = defaultdict(list)
    with open("./data/antibacterials.csv", "r") as fo:
        for line in tqdm(fo):
            _, group_id, smiles = line.strip().split(",")
            data[group_id].append(mol_to_barcode(smiles_to_mol(smiles), num_bits, radius)[0])
    
    # Split data
    antibacterial_train, antibacterial_test = train_test_split(data["antibacterial"], test_size=test_size, random_state=random_state, shuffle=True)
    not_antibacterial_train, not_antibacterial_test = train_test_split(data["not_antibacterial"], test_size=test_size, random_state=random_state, shuffle=True)
    
    X_train = np.array(antibacterial_train + not_antibacterial_train)
    y_train = np.array([1 for _ in antibacterial_train] + [0 for _ in not_antibacterial_train])
    X_test = np.array(antibacterial_test + not_antibacterial_test)
    y_test = np.array([1 for _ in antibacterial_test] + [0 for _ in not_antibacterial_test])

    # Create mask
    antibacterial_barcode = np.sum(antibacterial_train, axis=0)
    not_antibacterial_barcode = np.sum(not_antibacterial_train, axis=0)

    _, significant_bits = bit_enrichment(
        np.array((antibacterial_barcode, not_antibacterial_barcode)), 
        [len(antibacterial_train), len(not_antibacterial_train)],
        alpha,
        mtc=mtc)
    
    mask = np.vectorize(lambda x: x != 0)(np.sum(significant_bits, axis=0))
    print(f"Found `{sum(mask)}` significant features")

    # Apply mask to training and test sets
    X_train_masked = X_train[:, mask]
    X_test_masked = X_test[:, mask]

    # Train masked model
    model = RandomForestClassifier(n_estimators=1000, max_features="sqrt", random_state=random_state)
    model.fit(X_train_masked, y_train)
    test_preds = model.predict(X_test_masked)
    target_names = ["not_antibacterial_masked", "antibacterial_masked"]
    print(classification_report(y_test, test_preds, target_names=target_names))

    # Train unmasked model as baseline comparison
    model = RandomForestClassifier(n_estimators=1000, max_features="sqrt", random_state=random_state)
    model.fit(X_train, y_train)
    test_preds = model.predict(X_test)
    target_names = ["not_antibacterial", "antibacterial"]
    print(classification_report(y_test, test_preds, target_names=target_names))

    exit(0)


if __name__ == "__main__":
    main()
