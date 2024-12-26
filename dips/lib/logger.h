#include <ctime>
class Logger {
   public:
    static int mode;  // 1 for debug, 0 for info
    Logger() {}
    ~Logger();

    static void getDatetime() {
        time_t now = time(0);
        std::cout << ctime(&now) << std::endl;
    }
    template <typename _Ty, typename _Tx>
    static inline void LogInfo(_Tx key, _Ty val) {
        std::cout << key << ": " << val << std::endl;
    }

    template <typename _Ty>
    static inline void LogInfo(_Ty val) {
        std::cout << val << std::endl;
    }

    template <typename _Ty, typename _Tx>
    static inline void LogDebug(_Tx key, _Ty val) {
        if (mode == 1) std::cout << key << ": " << val << std::endl;
    }

    template <typename _Ty>
    static inline void LogError(_Ty val) {
        std::cout << val << std::endl;
    }

    template <typename _Ty, typename _Tx>
    static inline void LogError(_Tx key, _Ty val) {
        std::cout << key << ": " << val << std::endl;
    }

    template <typename _Ty>
    static inline void LogDebug(_Ty val) {
        if (mode == 1) std::cout << val << std::endl;
    }
};
