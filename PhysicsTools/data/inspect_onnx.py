import onnx
from onnx import shape_inference

model_path = "/work/gcelotto/btv_new/CMSSW_15_0_6/src/PhysicsTools/data/submod_out128_hyper_1802.onnx"

# Load model
model = onnx.load(model_path)

# Try shape inference (helps recover tensor shapes)
model = shape_inference.infer_shapes(model)

graph = model.graph

print("\n===== MODEL SUMMARY =====")

# Inputs
print("\n-- Inputs --")
for inp in graph.input:
    print(f"Name: {inp.name}")
    print(f"Type: {inp.type}")
    print()

# Outputs
print("\n-- Outputs --")
for out in graph.output:
    print(f"Name: {out.name}")
    print(f"Type: {out.type}")
    print()

# Nodes (core computation graph)
print("\n-- Nodes --")
for i, node in enumerate(graph.node):
    print(f"{i:04d}: {node.op_type}")
    print(f"  Inputs : {list(node.input)}")
    print(f"  Outputs: {list(node.output)}")
    print()

# Initializers (weights)
print("\n-- Initializers (weights/biases) --")
for init in graph.initializer:
    dims = list(init.dims)
    print(f"{init.name} | shape={dims} | dtype={init.data_type}")

print("\n===== END =====")