# Multimodal Deep Learning Model for T-Cell Classification

## Overview

This repository contains a Jupyter Notebook (`Multimodal-deep-learning_Model.ipynb`) implementing a multimodal deep learning model for classifying T-cells as "Neo" (neoantigen-specific) or "nonNeo" based on single-cell RNA sequencing data. The model integrates multiple data modalities, including:

- Gene expression profiles (18,314 genes).
- TCR (T-cell receptor) sequences (CDR3 alpha and beta chains).
- V(D)J gene usage (TRAV, TRAJ, TRBV, TRBD, TRBJ).
- Cell type annotations.
- HPV infection status.

The model uses embeddings for categorical features, Transformer encoders for sequence data, and fully connected layers for gene expression, fused via multi-head attention for binary classification.

Key features:
- Patient-level train-test splitting to avoid data leakage.
- 5-fold cross-validation with metrics (ROC AUC, PR AUPRC, Accuracy, F1-score).
- Prediction on unknown/unlabeled cells.

This model was developed for analyzing T-cell data in the context of HPV-related cancers, achieving an average ROC AUC of 0.9920 and F1-score of 0.9247 across folds.

## Prerequisites

### Environment
- Python 3.9 (or compatible version).
- Recommended: Create a Conda environment named `deep_learning_TR`:
  ```
  conda create -n deep_learning_TR python=3.9
  conda activate deep_learning_TR
  ```

### Dependencies
Install the required libraries using pip:
```
pip install pandas torch numpy matplotlib scikit-learn seaborn tqdm
```

Additional libraries from the code (pre-installed in many environments):
- `torch.utils.data` (for DataLoader and Dataset).
- `sklearn.model_selection` (for splitting).
- `sklearn.metrics` (for evaluation).

No internet access is needed during execution, as the code uses pre-loaded data and avoids external package installations.

## Data Preparation

### Input Data (located at `https://zenodo.org/uploads/17217824 `)
- **Main Dataset**: `Model_construction_data.csv`.
  - Columns: `ID`, `Label` (Neo/nonNeo), `cdr3_aa1` (alpha chain), `cdr3_aa2` (beta chain), `CTgene` (V(D)J genes), `HPV` (infection status), `Celltype`, `pMT`, and 18,314 gene expression columns.
  - Rows: 44,692 samples.
- **Metadata**: `input_meta-Molecular_cell.csv` (for patient IDs).
- **Unknown Data**: `unknown_cells_input.csv` (for predictions; similar format without labels).

### Preprocessing
- The notebook handles loading, merging patient metadata, and patient-level splitting.
- Features are encoded: categorical (one-hot/embeddings), sequences (amino acid properties), gene expressions (normalized floats).
- Outputs: Train/test indices, encoded tensors.

## Model Architecture

The `TCellClassifier` is a PyTorch-based neural network:

- **Embeddings**:
  - TCR genes (TRAV, TRAJ, TRBV, TRBD, TRBJ): Embedded into 32 dimensions.
  - Cell type: Embedded into 8 dimensions.
  - HPV: Embedded into 2 dimensions.
  - Amino acid sequences: Properties embedded into 32 dimensions.

- **Transformers**:
  - Two Transformer encoders (2 layers each) for alpha and beta CDR3 sequences.

- **Gene Expression Branch**:
  - Sequential FC layers: 18314 → 512 → 256 → 128 → 128 (with BatchNorm, ReLU, Dropout).

- **Fusion**:
  - Concatenate all features and apply multi-head attention.
  - FC layers: 364 → 128 → 64 → 2 (softmax for binary classification).

- **Loss and Optimizer**: Cross-Entropy Loss, Adam optimizer (lr=1e-5, weight_decay=1e-5).
- **Training**: 50 epochs, batch size 128, dropout rates 0.3-0.5.

## Usage

### Running the Notebook
1. Clone the repository:
   ```
   git clone https://github.com/your-username/multimodal-deep-learning-model.git
   cd multimodal-deep-learning-model
   ```
2. Adjust file paths in the notebook (e.g., CSV locations).

3. Run the notebook:
   - Use Jupyter: `jupyter notebook Multimodal-deep-learning_Model.ipynb`.
   - Or VS Code/Jupyter Lab.

### Key Steps in the Notebook
1. **Data Loading and Splitting**:
   - Load CSV, merge patient metadata.
   - Patient-level random split (80/20) with balanced positive ratios.

2. **Preprocessing**:
   - Function `pretreatment()`: Encodes features into tensors.

3. **Dataset and Dataloader**:
   - Custom `TCellDataset` and DataLoader.

4. **Cross-Validation**:
   - 5 folds, trains models, evaluates metrics.
   - Saves best model per fold (`best_model_fold_X.pth`).

5. **Prediction**:
   - Loads best model (based on F1-score).
   - Predicts on `unknown_cells_input.csv`.
   - Outputs: `predicted_unknown_cells.csv` (ID and predicted_label: Neo/nonNeo).

### Example Output
- Cross-validation summary:
  ```
  Average ROC AUC: 0.9920
  Average PR AUPRC: 0.9784
  Average test accuracy: 0.9388
  Average test F1: 0.9247
  ```
- Best fold: Fold 2 (F1: 0.9764).

## Evaluation Metrics
- **ROC AUC**: Measures discrimination.
- **PR AUPRC**: Precision-Recall AUC for imbalanced data.
- **Accuracy and F1-Score**: Overall performance.
- Visualizations (commented): ROC/PR curves, confusion matrices, loss plots (saved as SVG).

## Limitations
- Assumes balanced patient splitting; adjust `test_size` and `max_iter` for larger datasets.
- Hardcoded paths; modify for portability.
- No hyperparameter tuning; experiment with lr, dropout, etc.
- GPU recommended (uses CUDA if available).

## Contributing
Feel free to open issues or pull requests for improvements, such as adding hyperparameter tuning or supporting more modalities.

## License
This project is licensed under the MIT License.

## Contact
For questions, contact [wenzhishi@sina.com] or open an issue on GitHub.
