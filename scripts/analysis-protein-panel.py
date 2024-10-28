import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, roc_auc_score
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, confusion_matrix
from matplotlib import pyplot as plt


def main():
    # Data Preprocessing
    biomarker_clean = pd.read_csv('../data/biomarker_clean.csv')  # read data
    biomarker_clean = biomarker_clean.drop(
        ['Unnamed: 0'], axis=1)  # drop the unnamed column
    biomarker_clean[biomarker_clean['ados'].isna(
    )]['group'].value_counts()  # drop the ados column
    biomarker_clean['group'] = biomarker_clean['group'].apply(  # transfrom ASD to 1 and TD to 0
        lambda x: 1 if x == 'ASD' else 0).astype(int)

    # Initial Splitting
    X = biomarker_clean.loc[:, 'CHIP':]
    y = biomarker_clean.loc[:, 'group']

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, stratify=y, random_state=198)

    # Correlation
    corr_matrix = X_train.corr()
    to_drop = set()
    high_corr_pairs = [
        (col1, col2) for col1 in corr_matrix.columns for col2 in corr_matrix.columns if col1 != col2 and abs(corr_matrix.loc[col1, col2] > 0.9)
    ]

    for col1, col2 in high_corr_pairs:
        if col1 not in to_drop and col2 not in to_drop:
            to_drop.add(col2)

    # drop one of the correlated predictor
    X_train = X_train.drop(columns=to_drop)
    # drop one of the correlated predictor
    X_test = X_test.drop(columns=to_drop)

    # Recursive Feature Elimination in Logistic Regression
    log_reg_forward_model = LogisticRegression(
        solver='liblinear', random_state=197)
    rfe = RFE(estimator=log_reg_forward_model,
              n_features_to_select=20).fit(X_train, y_train)
    selected_features = X_train.columns[rfe.support_]

    # Random Forest (5-folds cross validation)
    rf_model = RandomForestClassifier(n_estimators=100,
                                      criterion='gini',
                                      max_depth=15,
                                      min_samples_split=5,
                                      bootstrap=True,
                                      random_state=197)

    # 5-folds cross validation
    kf = StratifiedKFold(n_splits=5, shuffle=True, random_state=197)
    feature_importance_rf = np.zeros(X_train.shape[1])

    for train_idx, test_idx in kf.split(X_train, y_train):
        X_kf_train, X_kf_test = X_train.iloc[train_idx], X_train.iloc[test_idx]
        y_kf_train, y_kf_test = y_train.iloc[train_idx], y_train.iloc[test_idx]

        rf_model.fit(X_kf_train, y_kf_train)

        feature_importance_rf += rf_model.feature_importances_

        # compare predicted ASD probability with outcome
        y_pred_prob = rf_model.predict_proba(X_kf_test)[:, 1]
        auc = roc_auc_score(y_kf_test, y_pred_prob)  # rocauc for each fold

    feature_importance_rf /= kf.get_n_splits()  # compute average feature importance

    # Select Top 20 Bio-markers Candidates
    feature_importance_df = pd.DataFrame({
        'Features': X_train.columns,
        'Importance': feature_importance_rf
    })
    top_20_features_rf = feature_importance_df.sort_values(
        by='Importance', ascending=False).head(20)

    overlapping_features = set(selected_features).union(  # union all selected features
        top_20_features_rf.loc[:, 'Features'])
    overlapping_features = list(overlapping_features)

    # Recursive Feature Elimination in Logistic Regression
    X_final_training = X_train.loc[:, overlapping_features]

    lg_model_final = LogisticRegression(solver='saga', penalty='elasticnet',
                                        C=1.0, l1_ratio=0.5, random_state=197, max_iter=500).fit(X_final_training, y_train)
    rfe = RFE(estimator=lg_model_final, n_features_to_select=6).fit(
        X_final_training, y_train)
    selected_features = X_final_training.columns[rfe.support_]
    print("Final Features:", selected_features)

    # Final model construction
    lg_model_final = LogisticRegression(solver='liblinear', random_state=197, max_iter=500).fit(
        X_train.loc[:, selected_features], y_train)

    # Plot Confusion Matrix
    cm = confusion_matrix(y_test, lg_model_final.predict(
        X_test.loc[:, selected_features]))
    fig, ax = plt.subplots(figsize=[8, 8])
    ax.imshow(cm, cmap="Blues")
    ax.xaxis.set(ticks=(0, 1), ticklabels=('Predicted TD', 'Predicted ASD'))
    ax.yaxis.set(ticks=(0, 1), ticklabels=('Actual TD', 'Actual ASD'))
    for i in range(2):
        for j in range(2):
            ax.text(j, i, cm[i, j], ha='center', va='center',
                    color='black', fontdict={'fontsize': 20})
    ax.grid(False)
    ax.set_title('Confusion Matrix of Final Predictions')
    plt.show()

    # Print evaluation report
    print(classification_report(y_test, lg_model_final.predict(
        X_test.loc[:, selected_features])))
    # Print ROCAUC score
    roc_auc_score(y_test, lg_model_final.predict(
        X_test.loc[:, selected_features]))


if __name__ == "__main__":
    main()
