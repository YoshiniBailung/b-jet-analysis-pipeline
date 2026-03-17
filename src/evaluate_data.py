def evaluate_data(model, df, config):

    features_for_train = config["features"]

    X = df[features_for_train]

    scores = model.predict(X, output_margin=False)

    return scores