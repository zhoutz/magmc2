#include <cstdlib>
#include <print>
#include <source_location>
#include <string_view>

inline
void
CHECK(bool condition, std::string_view msg,
      std::source_location const loc = std::source_location::current()) {
  if (condition)
    return;
  std::println(stderr, "[{}:{}] CHECK failed: {}", loc.file_name(), loc.line(),
               msg);
  std::abort();
}
