#pragma once

#include <cuda_runtime_api.h>

#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>

namespace medial_axis_3d::cuda_detail {

inline void check_cuda(cudaError_t error, const char* operation) {
    if (error == cudaSuccess) {
        return;
    }
    throw std::runtime_error(
        std::string(operation) + ": " + cudaGetErrorString(error)
    );
}

template <typename T>
class DeviceBuffer {
public:
    DeviceBuffer() = default;

    ~DeviceBuffer() {
        if (data_ != nullptr) {
            cudaFree(data_);
        }
    }

    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;

    DeviceBuffer(DeviceBuffer&& other) noexcept
        : data_(std::exchange(other.data_, nullptr)),
          capacity_(std::exchange(other.capacity_, 0)) {}

    DeviceBuffer& operator=(DeviceBuffer&& other) noexcept {
        if (this == &other) {
            return *this;
        }
        if (data_ != nullptr) {
            cudaFree(data_);
        }
        data_ = std::exchange(other.data_, nullptr);
        capacity_ = std::exchange(other.capacity_, 0);
        return *this;
    }

    void reserve(std::size_t count) {
        if (count <= capacity_) {
            return;
        }

        T* replacement = nullptr;
        check_cuda(
            cudaMalloc(
                reinterpret_cast<void**>(&replacement),
                count * sizeof(T)
            ),
            "cudaMalloc"
        );
        if (data_ != nullptr) {
            const cudaError_t free_error = cudaFree(data_);
            if (free_error != cudaSuccess) {
                cudaFree(replacement);
                check_cuda(free_error, "cudaFree");
            }
        }
        data_ = replacement;
        capacity_ = count;
    }

    T* data() noexcept {
        return data_;
    }

    const T* data() const noexcept {
        return data_;
    }

    std::size_t capacity() const noexcept {
        return capacity_;
    }

private:
    T* data_{nullptr};
    std::size_t capacity_{0};
};

}  // namespace medial_axis_3d::cuda_detail
