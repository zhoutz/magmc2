#include <cstdlib>
#include <print>
#include <source_location>
#include <string_view>

inline void
check(bool condition, std::string_view msg,
      std::source_location const loc = std::source_location::current()) {
  if (!condition) {
    std::println(stderr, "[{}:{}] check failed: {}", loc.file_name(),
                 loc.line(), msg);
    std::abort();
  }
}
