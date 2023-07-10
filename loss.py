import tensorflow as tf

def discriminative_loss_single(correct_label,prediction, feature_dim=3, label_shape=[512,256], 
                            delta_v=0.5, delta_d=3, param_var=1, param_dist=1, param_reg=0.001):
    
    ''' Discriminative loss for a single prediction/label pair.
    :param prediction: inference of network
    :param correct_label: instance label
    :feature_dim: feature dimension of prediction
    :param label_shape: shape of label
    :param delta_v: cutoff variance distance
    :param delta_d: curoff cluster distance
    :param param_var: weight for intra cluster variance
    :param param_dist: weight for inter cluster distances
    :param param_reg: weight regularization
    '''

    ### Reshape so pixels are aligned along a vector
    correct_label = tf.reshape(correct_label, [label_shape[1]*label_shape[0]])
    reshaped_pred = tf.reshape(prediction, [label_shape[1]*label_shape[0], feature_dim])

    ### Count instances
    unique_labels, unique_id, counts = tf.unique_with_counts(correct_label)
    counts = tf.cast(counts, tf.float32)
    num_instances = tf.size(unique_labels)

    segmented_sum = tf.math.unsorted_segment_sum(reshaped_pred, unique_id, num_instances)

    mu = tf.divide(segmented_sum, tf.reshape(counts, (-1, 1))) ### change tf.div=>tf.divide
    mu_expand = tf.gather(mu, unique_id)

    ### Calculate l_var
    distance = tf.norm(tf.subtract(mu_expand, reshaped_pred), axis=1)
    distance = tf.subtract(distance, delta_v)
    distance = tf.clip_by_value(distance, 0., distance)
    distance = tf.square(distance)

    l_var = tf.math.unsorted_segment_sum(distance, unique_id, num_instances)
    l_var = tf.divide(l_var, counts) ### change tf.div=>tf.divide
    l_var = tf.reduce_sum(l_var)
    l_var = tf.divide(l_var, tf.cast(num_instances, tf.float32))
    
    ### Calculate l_dist
    
    # Get distance for each pair of clusters like this:
    #   mu_1 - mu_1
    #   mu_2 - mu_1
    #   mu_3 - mu_1
    #   mu_1 - mu_2
    #   mu_2 - mu_2
    #   mu_3 - mu_2
    #   mu_1 - mu_3
    #   mu_2 - mu_3
    #   mu_3 - mu_3

    mu_interleaved_rep = tf.tile(mu, [num_instances, 1])
    mu_band_rep = tf.tile(mu, [1, num_instances])
    mu_band_rep = tf.reshape(mu_band_rep, (num_instances*num_instances, feature_dim))

    mu_diff = tf.subtract(mu_band_rep, mu_interleaved_rep)
    
    # Filter out zeros from same cluster subtraction
    intermediate_tensor = tf.reduce_sum(tf.abs(mu_diff),axis=1)
    zero_vector = tf.zeros(1, dtype=tf.float32)
    bool_mask = tf.not_equal(intermediate_tensor, zero_vector)
    mu_diff_bool = tf.boolean_mask(mu_diff, bool_mask)

    mu_norm = tf.norm(mu_diff_bool, axis=1)
    mu_norm = tf.subtract(2.*delta_d, mu_norm)
    mu_norm = tf.clip_by_value(mu_norm, 0., mu_norm)
    mu_norm = tf.square(mu_norm)

    l_dist = tf.reduce_mean(mu_norm)

    ### Calculate l_reg
    l_reg = tf.reduce_mean(tf.norm(mu, axis=1))

    param_scale = 1.
    l_var = param_var * l_var
    l_dist = param_dist * l_dist
    l_reg = param_reg * l_reg

    loss = param_scale*(l_var + l_dist + l_reg)
    
    return loss#, l_var, l_dist, l_reg
def discriminative_loss(correct_label,prediction, feature_dim=3, image_shape=[512,256], 
                            delta_v=0.5, delta_d=3, param_var=1, param_dist=1, param_reg=0.001):
    ''' Iterate over a batch of prediction/label and cumulate loss
    :return: discriminative loss
    '''
    def cond(label, batch, out_loss, i):
        return tf.less(i, tf.shape(batch)[0])
    def body(label, batch, out_loss, i):
        disc_loss = discriminative_loss_single(correct_label[i],prediction[i], feature_dim, image_shape, 
                        delta_v, delta_d, param_var, param_dist, param_reg)
        out_loss = out_loss.write(i, disc_loss)
        return label, batch, out_loss, i + 1
    # TensorArray is a data structure that support dynamic writing
    output_ta_loss = tf.TensorArray(dtype=tf.float32,
                   size=0,
                   dynamic_size=True)

    _, _, out_loss_op, _  = tf.while_loop(cond, body, [correct_label, 
                                                        prediction, 
                                                        output_ta_loss,
                                                        0])
    out_loss_op = out_loss_op.stack()
    disc_loss = tf.reduce_mean(out_loss_op)
    return disc_loss
def segmentation_loss_single(correct_label,prediction, eps=1.02,threshold=1e-10):
    correct_label = tf.reshape(correct_label, (-1,))
    correct_label=tf.cast(correct_label,tf.float32)
    reshaped_pred = tf.reshape(prediction, (-1,))  
    reshaped_pred=tf.cast(reshaped_pred,tf.float32)
    N=tf.size(correct_label)
    N=tf.cast(N,tf.float32)
    N_0=tf.size(tf.where(correct_label==0))
    N_0=tf.cast(N_0,tf.float32)
    N_1=tf.subtract(N,N_0)
    w_0=tf.divide(1.,tf.math.log(eps+tf.divide(N_0,N)))
    w_0=tf.cast(w_0, tf.float32)
    w_1=tf.divide(1.,tf.math.log(eps+tf.divide(N_1,N)))
    w_1=tf.cast(w_1, tf.float32)
    pred_threshold=tf.clip_by_value(reshaped_pred, threshold, reshaped_pred)
    one_pred_threshold=tf.clip_by_value(1-reshaped_pred, threshold, 1-reshaped_pred)
    loss_vector=tf.multiply(w_1,tf.multiply(correct_label,tf.math.log(pred_threshold)))\
        +tf.multiply(w_0,tf.multiply(1-correct_label,tf.math.log(one_pred_threshold)))
    loss=-tf.reduce_mean(loss_vector)
    #loss=tf.cast(loss,tf.float32)
    return loss
def segmentation_loss(correct_label,prediction, feature_dim=3, eps=1.02):
    ''' Iterate over a batch of prediction/label and cumulate loss
    :return: segmentation loss
    '''
    def cond(label, batch, out_loss, i):
        return tf.less(i, tf.shape(batch)[0])

    def body(label, batch, out_loss, i):
        seg_loss = segmentation_loss_single(correct_label[i],prediction[i])

        out_loss = out_loss.write(i, seg_loss)
        return label, batch, out_loss, i + 1
    # TensorArray is a data structure that support dynamic writing
    output_ta_loss = tf.TensorArray(dtype=tf.float32,
                   size=0,
                   dynamic_size=True)
    _, _, out_loss_op, _  = tf.while_loop(cond, body, [correct_label, 
                                                        prediction, 
                                                        output_ta_loss, 
                                                        0])
    out_loss_op = out_loss_op.stack()
    seg_loss = tf.reduce_mean(out_loss_op)
    return seg_loss