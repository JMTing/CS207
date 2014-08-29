#include <cstdio>
#include <functional>   // std::greater<T>
#include <vector>       // std::vector<T>

/** Current class implementation of binary_search
 *
 * TODO: Specification, generalization, and simplification
 */
template <typename T>
int binary_search(T* a, int low, int high, T v) {
  while (low <= high) {
    int mid = low + (high - low) / 2;
    if (a[mid] < v)
      low = mid + 1;
    else if (v < a[mid])
      high = mid - 1;
    else
      return mid;
  }
  return -1;
}

/** Current class implementation of binary_search
 *
 * TODO: Specification, generalization, and simplification
 */
template <typename T, typename Comparator>
int binary_search(T* a, int low, int high, T v, Comparator comp) {
  while (low <= high) {
    int mid = low + (high - low) / 2;
    if (comp(a[mid], v))
      low = mid + 1;
    else if (comp(v, a[mid]))
      high = mid - 1;
    else
      return mid;
  }
  return -1;
}


/** A simple function that returns @a a > @a b
 * Usage:
 * bool result = f_greater(4, 7);
 */
bool f_greater(int a, int b) {
  return a > b;
}

/** A function object (or functor) that can be called to return @a a > @a b
 * Usage:
 * s_greater gtr;   // Declare and construct an object
 * bool result = gtr(4, 7);
 *
 * Note that gtr(4, 7) is equivalent to gtr.operator()(4, 7)
 * This is just like (a < b) is equivalent to operator<(a, b)
 **/
struct s_greater {
  bool operator()(int a, int b) const {
    return a > b;
  }
};

int main()
{
  // Initialize a vector sorted in increasing order
  std::vector<int> v1 = {8, 18, 28, 36, 46, 65, 73, 83, 91};

  // Default implementation uses operator< and only works with ascending arrays
  int idx = binary_search(v1.data(), 0, v1.size(), 18);

  // Abstracted implementation uses a comparison function rather than op<
  // binary_search will work with descending arrays if we use > as our comparison

  // Initialize a vector sorted in decreasing order
  std::vector<int> v2 = {91, 83, 73, 65, 56, 36, 28, 18, 8};

  // Below are four (nearly) equivalent forms of passing a callable function
  // that returns a > b for two ints, a and b.

  // 1. Plain function parameter
  int f_idx = binary_search(v2.data(), 0, v2.size(), 18, f_greater);
  // 2. Function object (functor) parameter -- must construct an instance
  int s_idx = binary_search(v2.data(), 0, v2.size(), 18, s_greater());
  // 3. Standard library function object parameter
  int g_idx = binary_search(v2.data(), 0, v2.size(), 18, std::greater<int>());
  // 4. Lambda (unnamed) function parameter
  int l_idx = binary_search(v2.data(), 0, v2.size(), 18,
                            [](int a, int b) { return a > b; });

  // Print results
  printf("v1[%d] = %d\n", idx,   v1[idx]);
  printf("v2[%d] = %d\n", f_idx, v2[f_idx]);
  printf("v2[%d] = %d\n", s_idx, v2[s_idx]);
  printf("v2[%d] = %d\n", g_idx, v2[g_idx]);
  printf("v2[%d] = %d\n", l_idx, v2[l_idx]);

  // Fancy stuff!

  // v2 is also "sorted" by the last digit of the elements (increasing)
  // Let's find one that ends in 5
  int k1 = binary_search(v2.data(), 0, v2.size(), 05,
                         [](int a, int b) { return (a % 10) < (b % 10); });

  // v2 is also sorted by the first digit of the elements (decreasing)
  // Let's find one that begins in 5
  int k2 = binary_search(v2.data(), 0, v2.size(), 50,
                         [](int a, int b) { return (a / 10) > (b / 10); });

  printf("v2[%d] = %d\n", k1, v2[k1]);
  printf("v2[%d] = %d\n", k2, v2[k2]);

  // Note that this modularity comes at no cost.
  // We never modified or copied the data!
  // All operations are O(log(N)) and are as fast as if we had written
  //   a custom search function to accomplish each task above.

  // Wow, binary_search is already pretty powerful... Can it be improved?
  // Hint: Yes.
}
