###R Code for Fourth Year Thesis
####Prediction of Geomagnetically Induced Currents

Designing an early warning system for geomagnetically induced currents in power lines. Currently the data being used is simulated data.

Main file is <b>tf_test.R</b>, rest of files are functions to compute the spectral estimate of the data using the multitaper method, and then to compute the transfer function.

<b>get_spec_mtm.R</b> computes the spectral estimate of the data using the multitaper method. <b>get_spec.R</b> is an older version of the function that does not support blocked data.

<b>get_tf_all.R</b> computes the transfer function using the spectral estimates of the training data. <b>get_tf.R</b> is an older version of the function that treats each geomagnetic component individually, rather than as a whole.

<b>recon_func.R</b> is a function that reconstructs the frequency domain prediction back to time domain data. <b>reconstruction.R</b> was used to develop the recconstruction method. It does not support blocked data.