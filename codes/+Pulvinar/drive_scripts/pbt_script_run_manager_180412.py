


import sys
import os
import shutil
from server import Server
from lfads_wrapper.lfads_wrapper import lfadsWrapper
from lfads_wrapper.run_posterior_mean_sampling import run_posterior_sample_and_average
#import numpy as np

# parsing parameters from lfads-run-manager script
name = sys.argv[1]
run_save_path = sys.argv[2]     # where PBT will store the runs
data_dir = sys.argv[3]      # lfads data folder
lfads_output_dir = sys.argv[4]  # lfads

# Set PBT machine configuration

computers = [
    {'id': 'cortex',
     'ip': 'cortex.bme.emory.edu',
     'max_processes': 8,
     'process_start_cmd': '/snel/home/fzhu23/bin/PBT_HP_opt/pbt_opt/run_lfads_client.sh',
     'wait_for_process_start':True,
     },
    {'id': 'striate',
    'ip': 'striate.bme.emory.edu',
    'max_processes': 8,
    'process_start_cmd': '/snel/home/fzhu23/bin/PBT_HP_opt/pbt_opt/run_lfads_client.sh',
    'wait_for_process_start':True,
    },
]

# put the pbt run in the pbt_run folder
run_save_path = os.path.join(run_save_path, 'pbt_run')


''' create server object '''

svr = Server(name, run_save_path, num_workers=32, steps_to_ready=20,
             max_generations=50, exploit_method='binarytournament', exploit_param=[],
             explore_method='resample', explore_param=0.8, mode='minimize', force_overwrite=True,
             server_ip='localhost', port=27017, server_log_path=run_save_path, num_no_best_to_stop=5)


# add machines to be used by PBT
svr.add_computers(computers)

''' Specify the searchable parameters using explorable=True '''

# example if you want to use resample and specify the probablity of resampling per hyperparam
resample_param = {'sample_mode':'rand',
                  'sample_prob':0.5}

#svr.add_hp('keep_prob', (0.1, 1.0), init_sample_mode='rand', explore_method='resample', explore_param=0.3, limit_explore=True, explorable=True)
svr.add_hp('keep_prob', [0.95], explorable=False)

svr.add_hp('learning_rate_init', (0.0001, 0.05), init_sample_mode=[0.01], explore_method='resample', explore_param=0.3,  explorable=True)


''' Search L2 penalties '''

svr.add_hp('l2_gen_scale', (0.1, 10000), init_sample_mode='logrand', explorable=True)
svr.add_hp('l2_ic_enc_scale', (0.1, 10000), init_sample_mode='logrand', explorable=True)

svr.add_hp('l2_con_scale', (0.1, 10000), init_sample_mode='logrand', explorable=True)
svr.add_hp('l2_ci_enc_scale', (0.1, 10000), init_sample_mode='logrand', explorable=True)

svr.add_hp('l2_gen_2_factors_scale', (0, 10000), init_sample_mode=[00.0], explorable=False)
svr.add_hp('l2_ci_enc_2_co_in', (0, 10000), init_sample_mode=[00.0], explorable=False)

# KL costs
#svr.add_hp('kl_co_weight', (0.1, 10), init_sample_mode='logrand', explorable=True)
#svr.add_hp('kl_ic_weight', (0.1, 10), init_sample_mode='logrand', explorable=True)
svr.add_hp('kl_co_weight', [0.2], explorable=False)
svr.add_hp('kl_ic_weight', [0.2], explorable=False)


''' Path related params '''
svr.add_hp('data_dir', [data_dir], explorable=False)
svr.add_hp('data_filename_stem', ["lfads"], explorable=False)


''' Other fixed params '''
svr.add_hp('keep_ratio', [0.9], explorable=False)
svr.add_hp('max_ckpt_to_keep', [1], explorable=False)
svr.add_hp('max_ckpt_to_keep_lve', [1], explorable=False)
svr.add_hp('csv_log', ["fitlog"], explorable=False)
svr.add_hp('output_filename_stem', [""], explorable=False)
svr.add_hp('checkpoint_pb_load_name', ["checkpoint"], explorable=False)
svr.add_hp('checkpoint_name', ["lfads_vae"], explorable=False)
svr.add_hp('device', ["gpu:0"], explorable=False)
svr.add_hp('ps_nexamples_to_process', [100000000], explorable=False)


# important parameters
#svr.add_hp('val_cost_for_pbt', ['heldout_samp'], explorable=False)
svr.add_hp('val_cost_for_pbt', ['heldout_trial'], explorable=False)

svr.add_hp('cv_rand_seed', [0], explorable=False)
svr.add_hp('cv_keep_ratio', [1.0], explorable=False)


svr.add_hp('factors_dim', [40], explorable=False)
svr.add_hp('in_factors_dim', [40])
svr.add_hp('do_train_readin', [False])

svr.add_hp('ic_dim', [64], explorable=False)
svr.add_hp('ic_enc_dim', [64], explorable=False)
svr.add_hp('gen_dim', [64], explorable=False)
svr.add_hp('batch_size', [160], explorable=False)
svr.add_hp('allow_gpu_growth', [True], explorable=False)
svr.add_hp('learning_rate_decay_factor', [1], explorable=False)
svr.add_hp('learning_rate_stop', [0.00000001], explorable=False)
svr.add_hp('learning_rate_n_to_compare', [float('inf')], explorable=False)
svr.add_hp('do_reset_learning_rate', [False], explorable=False)
svr.add_hp('con_ci_enc_in_dim', [40], explorable=False)
svr.add_hp('con_fac_in_dim', [40], explorable=False)
svr.add_hp('max_grad_norm', [200.0], explorable=False)
svr.add_hp('cell_clip_value', [5.0], explorable=False)
svr.add_hp('output_dist', ['poisson'], explorable=False)

svr.add_hp('co_dim', [15], explorable=False)
svr.add_hp('do_causal_controller', [False], explorable=False)
svr.add_hp('controller_input_lag', [1], explorable=False)
svr.add_hp('ci_enc_dim', [64], explorable=False)
svr.add_hp('con_dim', [64], explorable=False)

svr.add_hp('ic_prior_var', [0.1], explorable=False)
svr.add_hp('ic_post_var_min', [0.0001], explorable=False)
svr.add_hp('co_prior_var', [0.1], explorable=False)
svr.add_hp('co_post_var_min', [0.0001], explorable=False)
svr.add_hp('do_feed_factors_to_controller', [True], explorable=False)
svr.add_hp('feedback_factors_or_rates', ['factors'])
#svr.add_hp('l2_gen_scale', [500.0])
#svr.add_hp('l2_con_scale', [500.0])
svr.add_hp('kl_increase_steps', [900])
svr.add_hp('temporal_spike_jitter_width', [0])
svr.add_hp('inject_ext_input_to_gen', [False])


''' Start PBT '''
best_worker = svr.start_pbt()

''' get the path to the best worker '''
best_worker_path = best_worker.run_save_path

''' do posterior mean sample and average on the best model '''
# PM bach size
# pass None to use the training batch size (not recommended)
PM_batch_size = 200

# run posterior sample and average on the best worker
run_posterior_sample_and_average(best_worker_path, PM_batch_size)

''' copy the best worker to lfadsOutput for run-manager to load '''

src = best_worker_path
dest = lfads_output_dir

# create lfadsOutput directory if doesn't exist
if os.path.exists(dest):
    shutil.rmtree(dest)

os.makedirs(dest)

src_files = os.listdir(src)
for file_name in src_files:
    full_file_name = os.path.join(src, file_name)
    if (os.path.isfile(full_file_name)):
        shutil.copy(full_file_name, dest)




