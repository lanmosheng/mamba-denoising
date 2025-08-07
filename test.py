import torch

# 加载保存的张量
input_data = torch.load("debug_input.pt")
output_data = torch.load("debug_output.pt")
label_data = torch.load("debug_label.pt")

print("🔎 input_data shape:", input_data.shape)
print("🔎 output_data shape:", output_data.shape)
print("🔎 label_data shape:", label_data.shape)

# 检查是否存在异常值
print("✅ input_data - NaN:", torch.isnan(input_data).any().item(), 
      "Inf:", torch.isinf(input_data).any().item(), 
      "Max:", input_data.max().item(), 
      "Min:", input_data.min().item(), 
      "Mean:", input_data.mean().item())

print("✅ label_data - NaN:", torch.isnan(label_data).any().item(), 
      "Inf:", torch.isinf(label_data).any().item(), 
      "Max:", label_data.max().item(), 
      "Min:", label_data.min().item(), 
      "Mean:", label_data.mean().item())

print("❌ output_data - NaN:", torch.isnan(output_data).any().item(), 
      "Inf:", torch.isinf(output_data).any().item())
with torch.no_grad():
    print("Output stats:")
    print(" - Max:", output_data.max().item())
    print(" - Min:", output_data.min().item())
    print(" - Mean:", output_data.mean().item())
    print(" - Std :", output_data.std().item())

    # 哪些位置是 NaN 或 Inf？
    print(" - NaN indices:", torch.nonzero(torch.isnan(output_data)))
    print(" - Inf indices:", torch.nonzero(torch.isinf(output_data)))
