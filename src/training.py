import xgboost as xgb
from sklearn.model_selection import train_test_split
from hipe4ml.model_handler import ModelHandler
import yaml

def load_params(pt_bin):

    with open("config/xgbparams.yaml") as f:
        cfg = yaml.safe_load(f)

    return cfg["pt_bins"][pt_bin]

def train_model(df, pt_bin, config):

    features = config["features"]
    label = config["label"]
    X = df[features]
    y = df[label]

    params = load_params(pt_bin)

    weights = None
    if config["training"]["use_weights"]:
        print("Training with Event Weights")
        weights = df[config["weight"]]/min(df[config["weight"]])
        print(params)
        X_train, X_test, y_train, y_test, w_train, w_test = train_test_split(
            X, y, weights,
            test_size=config["training"]["test_size"],
            random_state=42
        )
        model = xgb.XGBClassifier(**params)

        model.fit(
            X_train,
            y_train,
            sample_weight = w_train
        )
    else:
        print("Training without Event Weights")
        print(params)
        X_train, X_test, y_train, y_test = train_test_split(
            X, y,
            test_size=config["training"]["test_size"],
            random_state=42
        )

        model = xgb.XGBClassifier(**params)

        model.fit(
            X_train,
            y_train
        )

    model_hdl = ModelHandler(model, features, model_params=None)

    train_test_data = [X_train, y_train, X_test, y_test]

    return model_hdl, train_test_data