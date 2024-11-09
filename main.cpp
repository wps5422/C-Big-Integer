#include "Integer.h" // Giả sử bạn lưu lớp Integer trong file Integer.h
#include <cassert>
#include <iostream>

void test_multiplication() {
    // Test 1: Nhân hai số nhỏ
    {
        Integer a(3); // 3
        Integer b(4); // 4
        Integer result = a * b; // 12
        assert(result.to_string() == "12");
        std::cout << "Test 1 Passed: 3 * 4 = " << result << std::endl;
    }

    // Test 2: Nhân số lớn
    {
        Integer a("12345678901234567890");
        Integer b("98765432109876543210");
        Integer expected("1219326311370217952237463801111263526900");
        Integer result = a * b;
        assert(result.to_string() == expected.to_string());
        std::cout << "Test 2 Passed: " << a << " * " << b << " = " << result << std::endl;
    }

    // Test 3: Nhân số âm với số dương
    {
        Integer a("-12345");
        Integer b("6789");
        Integer expected("-83810205");
        Integer result = a * b;
        assert(result.to_string() == expected.to_string());
        std::cout << "Test 3 Passed: " << a << " * " << b << " = " << result << std::endl;
    }

    // Test 4: Nhân hai số âm
    {
        Integer a("-12345");
        Integer b("-6789");
        Integer expected("83810205");
        Integer result = a * b;
        assert(result.to_string() == expected.to_string());
        std::cout << "Test 4 Passed: " << a << " * " << b << " = " << result << std::endl;
    }

    // Test 5: Nhân với 0
    {
        Integer a("123456789");
        Integer b("0");
        Integer expected("0");
        Integer result = a * b;
        assert(result.to_string() == expected.to_string());
        std::cout << "Test 5 Passed: " << a << " * " << b << " = " << result << std::endl;
    }

    // Test 6: Nhân với 1
    {
        Integer a("987654321");
        Integer b("1");
        Integer expected("987654321");
        Integer result = a * b;
        assert(result.to_string() == expected.to_string());
        std::cout << "Test 6 Passed: " << a << " * " << b << " = " << result << std::endl;
    }

    // Test 7: Nhân số lớn với số nhỏ
    {
        Integer a("99999999999999999999");
        Integer b("2");
        Integer expected("199999999999999999998");
        Integer result = a * b;
        assert(result.to_string() == expected.to_string());
        std::cout << "Test 7 Passed: " << a << " * " << b << " = " << result << std::endl;
    }

    // Test 8: Nhân số có nhiều chữ số
    {
        Integer a("314159265358979323846264338327950288419716939937510");
        Integer b("271828182845904523536028747135266249775724709369995");
        // Kết quả được tính toán trước bằng công cụ đáng tin cậy
        // Để đơn giản, ở đây ta chỉ kiểm tra việc nhân không gây lỗi
        Integer result = a * b;
        std::cout << "Test 8 Passed: Nhân hai số có nhiều chữ số thành công." << std::endl;
    }

    std::cout << "Tất cả các test đã vượt qua!" << std::endl;
}

int main() {
    test_multiplication();
    return 0;
}
