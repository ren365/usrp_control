// set up our masks, defining the pin numbers
#define AMP_GPIO_MASK   (1 << 6) | (1 << 10)
#define MAN_GPIO_MASK   (1 << 4)
#define ATR_MASKS       (AMP_GPIO_MASK | MAN_GPIO_MASK)
// set up our values for ATR control: 1 for ATR, 0 for manual
#define ATR_CONTROL     (AMP_GPIO_MASK)
// set up the GPIO directions: 1 for output, 0 for input
#define GPIO_DDR        (AMP_GPIO_MASK | MAN_GPIO_MASK)
// assume an existing multi_usrp object, called "usrp"
// now, let's do the basic ATR setup
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/utils/safe_main.hpp>
#include <iostream>

int UHD_SAFE_MAIN(int argc, char* argv[])
{

    uhd::device_addr_t my_usrp_ip;
    my_usrp_ip["addr0"] = "192.168.10.2";
    uhd::usrp::multi_usrp::sptr  my_usrp = uhd::usrp::multi_usrp::make(my_usrp_ip);
    std::cout << boost::format("Using Device: %s") % my_usrp->get_pp_string() << std::endl;

    std::cout << "Available GPIO banks: " << std::endl;
    auto banks = my_usrp->get_gpio_banks(0);
    for (auto& bank : banks) {
        std::cout << "* " << bank << std::endl;
    }

    my_usrp->set_gpio_attr("FP0", "CTRL", ATR_CONTROL, ATR_MASKS);
    my_usrp->set_gpio_attr("FP0", "DDR", GPIO_DDR, ATR_MASKS);
    // let's manually set GPIO4 high
    my_usrp->set_gpio_attr("FP0", "OUT", MAN_GPIO_MASK, MAN_GPIO_MASK);
    // finally, let's set up GPIO6 as we described above
    my_usrp->set_gpio_attr("FP0", "ATR_0X", AMP_GPIO_MASK, AMP_GPIO_MASK);
    my_usrp->set_gpio_attr("FP0", "ATR_RX", 0, AMP_GPIO_MASK);
    my_usrp->set_gpio_attr("FP0", "ATR_TX", 0, AMP_GPIO_MASK);//AMP_GPIO_MASK
    // usually, you would want to also make this pin go high when doing
    // full-duplex, but not in this example
    my_usrp->set_gpio_attr("FP0", "ATR_XX", 0, AMP_GPIO_MASK);


    return EXIT_SUCCESS;
}

