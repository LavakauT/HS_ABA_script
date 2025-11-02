import os
import sys

def generate_command(folder_path):
    # 讀取並排序檔案
    files = sorted(f for f in os.listdir(folder_path) if f.endswith(".txt"))
    
    # 組成完整路徑
    input_files = [os.path.join(folder_path, f) for f in files]
    
    # 移除前後字串，生成 group 名稱
    group_names = []
    for f in files:
        name = f
        name = name.replace("upset_up_", "").replace("upset_down_", "")
        name = name.replace("_intersection_genes.txt", "")
        name = name.replace("_genes.txt", "")
        name = name.replace("_intersection.txt", "")
        name = name.replace(".txt", "")
        group_names.append(name)
    
    # 組成指令
    command = (
        "python src/go_term_analysis.py \\\n"
        "--input-files \\\n"
        + " \\\n".join(f'"{path}"' for path in input_files)
        + " \\\n--group-names \\\n"
        + " \\\n".join(f'"{name}"' for name in group_names)
    )
    return command

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("用法: python generate_command.py <資料夾路徑>")
        sys.exit(1)
    
    folder_path = sys.argv[1]
    if not os.path.isdir(folder_path):
        print(f"錯誤：找不到資料夾 {folder_path}")
        sys.exit(1)
    
    print(generate_command(folder_path))
