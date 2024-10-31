import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
from sklearn.feature_selection import RFE
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import roc_auc_score, classification_report


def main():
    ### Data Preprocessing ###
    biomarker_clean = pd.read_csv('../data/biomarker_clean.csv')  # read data

    biomarker_clean = biomarker_clean.drop(
        ['Unnamed: 0'], axis=1)  # drop the unnamed column
    biomarker_clean[biomarker_clean['ados'].isna(
    )]['group'].value_counts()  # drop the ados column
    biomarker_clean['group'] = biomarker_clean['group'].apply(  # transfrom ASD to 1 and TD to 0
        lambda x: 1 if x == 'ASD' else 0).astype(int)

    X = biomarker_clean.loc[:, 'CHIP':]
    y = biomarker_clean.loc[:, 'group']
    X_train, X_test, y_train, y_test = train_test_split( # initial splitting
        X, y, test_size=0.2, stratify=y, random_state=198)



    ### Correlation Feature Selection ###
    corr_matrix = X_train.corr()
    to_drop = set()
    high_corr_pairs = [
        (col1, col2) for col1 in corr_matrix.columns for col2 in corr_matrix.columns if col1 != col2 and abs(corr_matrix.loc[col1, col2] > 0.9)
    ]

    for col1, col2 in high_corr_pairs:
        if col1 not in to_drop and col2 not in to_drop:
            to_drop.add(col2)

    X_train = X_train.drop(columns=to_drop) # drop from training set
    X_test = X_test.drop(columns=to_drop) # drop from testing set



    ### Recursive Feature Elimination in Logistic Regression - L1 ###
    log_reg_forward_model = LogisticRegression( # logistic regression initializer
        solver='liblinear', 
        random_state=197)
    
    rfe = RFE(estimator=log_reg_forward_model,
              n_features_to_select=20) # rfe initializer
    rfe.fit(X_train, y_train)

    selected_features = X_train.columns[rfe.support_]
    print("Selected features by step 2:", selected_features)



    ### Random Forest - 5 folds CV ###
    rf_model = RandomForestClassifier(n_estimators=100,
                                      criterion='gini',
                                      max_depth=15,
                                      min_samples_split=5,
                                      bootstrap=True,
                                      random_state=197)
    
    kf = StratifiedKFold(n_splits=5, 
                         shuffle=True, 
                         random_state=197)
    
    feature_importance_rf = np.zeros(X_train.shape[1])
    for train_idx, test_idx in kf.split(X_train, y_train):
        X_kf_train = X_train.iloc[train_idx]
        y_kf_train = y_train.iloc[train_idx]

        rf_model.fit(X_kf_train, y_kf_train)
        
        feature_importance_rf += rf_model.feature_importances_

    feature_importance_rf /= kf.get_n_splits()

    feature_importance_df = pd.DataFrame({
        'Features': X_train.columns,
        'Importance': feature_importance_rf
    })
    top_20_features_rf = feature_importance_df.sort_values(
        by='Importance', ascending=False).head(20)
    
    print("\nSelected features by step 3:", top_20_features_rf.columns)



    ### Find the Union features ###
    overlapping_features = set(selected_features).union(
        top_20_features_rf.loc[:, 'Features'])
    overlapping_features = list(overlapping_features)
    print("\nUnion features from step 2 and 3:", overlapping_features)



    ### Recursive Feature Elimination in Elastic Net Logistic Regression ###
    X_overlapping_features = X_train.loc[:, overlapping_features]

    lg_model_eln = LogisticRegression(solver='saga', 
                                      penalty='elasticnet',
                                      C=1.0, 
                                      l1_ratio=0.5, 
                                      random_state=197, 
                                      max_iter=500)
    
    sfs = RFE(estimator=lg_model_eln,
              n_features_to_select=6)
    sfs.fit(X_overlapping_features, y_train)

    selected_features = sfs.get_support(indices=True) # select features
    protein_panel_name = X_overlapping_features.columns[selected_features]
    X_protein_panel = X_overlapping_features.iloc[:, selected_features]

    print("\nFinal protein panel:", protein_panel_name)



    ### Final model construction ###
    lg_model_final = LogisticRegression(solver='liblinear',
                                        random_state=197,
                                        max_iter=500).fit(X_protein_panel, y_train)

    # Print evaluation report
    print(classification_report(y_test, lg_model_final.predict(
        X_test.loc[:, protein_panel_name])))

    # Print ROCAUC score
    print(roc_auc_score(y_test, lg_model_final.predict(
        X_test.loc[:, protein_panel_name])))


if __name__ == "__main__":
    main()
