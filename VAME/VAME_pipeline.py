import os
import shutil
import glob
import numpy as np
import subprocess
import vame

def run_vame_behavior_pipeline(project_name, working_dir, source_h5_dir):
    """Phase 1: VAME Preprocessing, Egocentric Alignment, and HMM Segmentation"""
    print(f"\n--- Starting Phase 1: VAME Behavioral Segmentation ---")
    
    # 1. Initialize Project & Move Pose Files
    h5_source_files = glob.glob(os.path.join(source_h5_dir, "**", "*.h5"), recursive=True)
    video_list = list(set([os.path.dirname(f) for f in h5_source_files]))
    
    if not video_list:
        print("No .h5 files found. Skipping VAME initialization.")
        return
        
    config_path = vame.init_new_project(
        project=project_name, 
        videos=video_list, 
        working_directory=working_dir, 
        videotype='.mp4'
    )
    
    dest_pose_dir = os.path.join(os.path.dirname(config_path), "videos", "pose_estimation")
    for h5_file in h5_source_files:
        shutil.copy(h5_file, dest_pose_dir)
        
    # 2. Egocentric Alignment
    print("Running Egocentric Alignment...")
    vame.egocentric_alignment(config_path)
    
    # 3. Crop Keypoints to Target Columns (snout, paws, tailbase)
    print("Cropping keypoints for Blackbox integration...")
    wanted_parts = ['snout', 'lfpaw', 'rfpaw', 'rhpaw', 'lhpaw', 'tailbase']
    all_parts = [
        'tailtip', 'tailbase', 'anal', 'genital', 'hip', 'sternumtail', 'sternumhead', 'neck', 'snout', 
        'lcheek', 'rcheek', 'lhip', 'rhip', 'lshoulder', 'rshoulder', 'lankle', 'lhpaw', 'lhpd1b', 
        'lhpd1t', 'lhpd2b', 'lhpd2t', 'lhpd3b', 'lhpd3t', 'lhpd4b', 'lhpd4t', 'lhpd5b', 'lhpd5t', 
        'lfpaw', 'lfpd1b', 'lfpd1t', 'lfpd2b', 'lfpd2t', 'lfpd3b', 'lfpd3t', 'lfpd4b', 'lfpd4t', 
        'rankle', 'rhpaw', 'rhpd1b', 'rhpd1t', 'rhpd2b', 'rhpd2t', 'rhpd3b', 'rhpd3t', 'rhpd4b', 
        'rhpd4t', 'rhpd5b', 'rhpd5t', 'rfpaw', 'rfpd1b', 'rfpd1t', 'rfpd2b', 'rfpd2t', 'rfpd3b', 
        'rfpd3t', 'rfpd4b', 'rfpd4t'
    ]
    
    indices_to_keep = []
    for part in wanted_parts:
        idx = all_parts.index(part)
        indices_to_keep.extend([idx*2, idx*2+1])
        
    data_root = os.path.join(os.path.dirname(config_path), 'data')
    for folder in os.listdir(data_root):
        folder_path = os.path.join(data_root, folder)
        if os.path.isdir(folder_path):
            npy_file = os.path.join(folder_path, f"{folder}-PE-seq.npy")
            if os.path.exists(npy_file):
                data = np.load(npy_file)
                if data.shape[0] == 114:
                    data = data.T
                cropped_data = data[:, indices_to_keep].T
                np.save(npy_file, cropped_data)
                
    # 4. Create Trainset & Train Model
    print("Training VAME Model...")
    vame.create_trainset(config_path)
    vame.train_model(config_path)
    
    # 5. Pose Segmentation (HMM)
    print("Executing HMM Pose Segmentation...")
    vame.pose_segmentation(config_path)
    print("Phase 1 Complete.")

def run_matlab_script(script_path):
    """Executes a MATLAB script headlessly from Python."""
    if not os.path.exists(script_path):
        print(f"Error: MATLAB script not found at {script_path}")
        return False
        
    print(f"\n--- Executing MATLAB Script: {os.path.basename(script_path)} ---")
    # -batch runs the script without opening the MATLAB GUI and exits automatically
    cmd = f'matlab -batch "run(\'{script_path}\')"'
    
    try:
        subprocess.run(cmd, shell=True, check=True)
        print(f"Successfully completed {os.path.basename(script_path)}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"MATLAB Execution Error for {os.path.basename(script_path)}: {e}")
        return False

if __name__ == "__main__":
    # =====================================================================
    # 1. CONFIGURATION PATHS
    # =====================================================================
    # VAME Paths
    PROJECT_NAME = 'SNI_Longitudinal_Full'
    WORKING_DIR = r"E:\Behavior\VAME_Projects"
    SOURCE_H5_DIR = r"E:\Behavior\Blackbox_analysis"
    
    # MATLAB Script Paths (Assumes they are saved in the same folder as this Python script)
    # Update these filenames to match how they are saved in your repository
    CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
    MATLAB_PHOTOMETRY_BATCH = os.path.join(CURRENT_DIR, "step1_photometry_batch.m")
    MATLAB_ALIGNMENT = os.path.join(CURRENT_DIR, "step2_vame_photometry_alignment.m")
    MATLAB_AUROC = os.path.join(CURRENT_DIR, "step3_auroc_permutation.m")
    MATLAB_STATS_RECOVERY = os.path.join(CURRENT_DIR, "step4_final_stats_recovery.m")

    # =====================================================================
    # 2. EXECUTE PIPELINE
    # =====================================================================
    print("Starting Automated Longitudinal Behavior & Neural Alignment Pipeline")
    
    # Run Phase 1: Python VAME Extraction
    run_vame_behavior_pipeline(PROJECT_NAME, WORKING_DIR, SOURCE_H5_DIR)
    
    # Run Phase 2: TDT Photometry Cleaning (MATLAB)
    # This will prompt you to select the parent directory for the TDT folders
    run_matlab_script(MATLAB_PHOTOMETRY_BATCH)
    
    # Run Phase 3: Synchronize VAME motifs with cleaned Z-scores (MATLAB)
    run_matlab_script(MATLAB_ALIGNMENT)
    
    # Run Phase 4: Permutation testing and Prism histogram exports (MATLAB)
    run_matlab_script(MATLAB_AUROC)
    
    # Run Phase 5: Hedges' g, KS tests, and manuscript summary generation (MATLAB)
    run_matlab_script(MATLAB_STATS_RECOVERY)
    
    print("\n======================================================")
    print("FULL PIPELINE EXECUTION COMPLETE.")
    print("Check the respective output directories for figures and Prism-ready Excel files.")
    print("======================================================")