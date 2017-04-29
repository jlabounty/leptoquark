from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import urllib
import urllib.request

import numpy as np
import tensorflow as tf

# Data sets
IRIS_TRAINING = "JetSummary_1000_training.csv"

IRIS_TEST = "JetSummary_1000_validation.csv"

def main():

  # Load datasets.
  training_set = tf.contrib.learn.datasets.base.load_csv_with_header(
      filename=IRIS_TRAINING,
      target_dtype=np.int,
      features_dtype=np.float32)
  test_set = tf.contrib.learn.datasets.base.load_csv_with_header(
      filename=IRIS_TEST,
      target_dtype=np.int,
      features_dtype=np.float32)

  # Specify that all features have real-value data
  feature_columns = [tf.contrib.layers.real_valued_column("", dimension=19)]

  # Build 3 layer DNN with 10, 20, 10 units respectively.
  classifier = tf.contrib.learn.DNNClassifier(feature_columns=feature_columns,
                                              hidden_units=[10, 20, 40, 20, 10],
                                              n_classes=2,
                                              model_dir="/tmp/jet_model")
  # Define the training inputs
  def get_train_inputs():
    x = tf.constant(training_set.data)
    y = tf.constant(training_set.target)

    return x, y

  # Fit model.
  classifier.fit(input_fn=get_train_inputs, steps=4000)

  # Define the test inputs
  def get_test_inputs():
    x = tf.constant(test_set.data)
    y = tf.constant(test_set.target)

    return x, y

  # Evaluate accuracy.
  accuracy_score = classifier.evaluate(input_fn=get_test_inputs,
                                       steps=1)["accuracy"]

  print("\nTest Accuracy: {0:f}\n".format(accuracy_score))

  # Classify two new flower samples.
  def new_samples():
    return np.array(
      [[280, 252, 196, 60, 10, 0, 0.215892, 0.203238, 0.739729, 0.183651, -0.0132682, 0.203238, 0.00341859, 0.183651, 0.000781234, 0.203723, 0.00141542, 0.183662, 55.2661],
       [57, 56, 52, 25, 8, 2, 2.82755, 0.279101, 1.86751, 0.22915, -0.143172, 0.279101, 0.0304653, 0.22915, -0.0571866, 0.292046, 0.0224339, 0.229291, 76.0819]], dtype=np.float32)

  predictions = list(classifier.predict(input_fn=new_samples))

  print(
      "New Samples, Class Predictions:    {}\n"
      .format(predictions))

if __name__ == "__main__":
    main()
