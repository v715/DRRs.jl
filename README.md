# DRRs

[![Build Status](https://github.com/v715/DRRs.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/v715/DRRs.jl/actions/workflows/CI.yml?query=branch%3Amain)

**NOTE: THIS PACKAGE IS NO LONGER MAINTAINED.**

*Parallelized, GPU-accelerated, and differentiable digitally reconstructed radiographs in Julia.*

⚠️ `DDRs.jl` was an experimental digitally reconstructed radiograph (DRR) generator in Julia. It's aim was to be a GPU-accelerated, auto-differentiable renderer (a la [RayTracer.jl](https://github.com/avik-pal/RayTracer.jl)) for medical images, but challenges using `CUDA.jl` and `Zygote.jl` made it difficult to reach a minimal viable prototype.

🫁 `DRRs.jl` inspired a PyTorch version that is fully functional and publically available at [`DiffDRR`](https://github.com/v715/DiffDRR)

🤷🏾‍♀️ Lessons learned from `DiffDRR` may be ported back to Julia at some point, but having to write custom differentiation rules is really annoying (compared to PyTorch's autograd, which just works 😕)
