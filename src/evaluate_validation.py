import numpy as np
from sklearn.metrics import roc_auc_score

def evaluate_validation(model, X_test, y_test):

    scores = model.predict(X_test, output_margin = False)

    roc_auc = roc_auc_score(y_test, scores, multi_class='ovo')

    print(f"ROC AUC Score: {roc_auc}")

    return scores