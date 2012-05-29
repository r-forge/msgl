/*
 * sgl_debug.h
 *
 *  Created on: Jul 19, 2011
 *      Author: martin
 */

#ifndef CONFIG_H_
#define CONFIG_H_

#ifdef SGL_CONVERGENCE_CHECK
#define CONVERGENCE_CHECK(limit) sgl::natural convergence_itr = 0; sgl::natural convergence_limit = limit;
#define CONVERGENCE_CHECK_INCREASE ++convergence_itr; if(convergence_itr > convergence_limit) throw std::runtime_error(create_error_msg("Convergence problems", __FILE__, __LINE__));
#else
#define CONVERGENCE_CHECK(limit)
#define CONVERGENCE_CHECK_INCREASE
#endif

#ifdef SGL_CATCH_EXCEPTIONS
#define SGL_EXCEPTIONS std::exception
#define SGL_EXCEPTIONS_GENRAL ...
#else
#define SGL_EXCEPTIONS dummy_exception_a
#define SGL_EXCEPTIONS_GENRAL dummy_exception_b
#endif

#ifndef SGL_WARNING
#define SGL_WARNING(msg) std::cerr << msg << std::endl;
#endif

#ifndef SGL_ERROR
#define SGL_ERROR(msg) std::cerr << msg << std::endl; exit(1);
#endif

#ifndef SGL_MSG
#define SGL_MSG(msg) std::cout << msg << std::endl;
#endif

#ifndef SGL_INTERRUPT_INIT
#define SGL_INTERRUPT_INIT
#endif

#ifndef SGL_INTERRUPT_CHECK
#define SGL_INTERRUPT_CHECK
#endif

#ifndef SGL_INTERRUPT_RESET
#define SGL_INTERRUPT_RESET
#endif

#ifndef SGL_INTERRUPT
#define SGL_INTERRUPT
#endif

//Debug

#ifdef SGL_RUNTIME_CHECKS
//TODO templates should be used for this
#define ASSERT_IS_FINITE(x) if(!sgl::is_finite(x)) throw std::runtime_error(create_error_msg("Numerical problems", __FILE__, __LINE__));
#define ASSERT_IS_NUMBER(x) if(x != x) throw std::runtime_error(create_error_msg("Numerical problems", __FILE__, __LINE__));
#define ASSERT_IS_POSITIVE(x) if(x <= 0) throw std::runtime_error(create_error_msg("Numerical problems", __FILE__, __LINE__));
#define ASSERT_IS_NON_NEGATIVE(x) if(x < 0) throw std::runtime_error(create_error_msg("Numerical problems", __FILE__, __LINE__));
#else
#define ASSERT_IS_FINITE(x) //do nothing
#define ASSERT_IS_NUMBER(x) //do nothing
#define ASSERT_IS_POSITIVE(x) //do nothing
#define ASSERT_IS_NON_NEGATIVE(x) //do nothing
#endif

#ifdef SGL_DEBUG_SIMPLE
//TODO templates should be used for this
#define ASSERT_IS_ZERO(x) if(accu(x != 0) > 0) throw std::runtime_error(create_error_msg("Numerical problems", __FILE__, __LINE__));
#define ASSERT_IS_NON_ZERO(x) if(accu(x != 0) == 0) throw std::runtime_error(create_error_msg("Numerical problems", __FILE__, __LINE__));
#else
#define ASSERT_IS_ZERO(x) //do nothing
#define ASSERT_IS_NON_ZERO(x) //do nothing
#endif

#ifdef SGL_DEBUG_COMPLEX
#define SGL_DEBUG_BLOCK_ACTIVE
#define SGL_DEBUG_GB
#endif

#ifdef SGL_DEBUG_INFO_ALL
#define SGL_DEBUG_INFO_GB_OPT
#define SGL_DEBUG_INFO_ACTIVE_SET
#define SGL_DEBUG_INFO_STEPSIZE
#define SGL_DEBUG_INFO_QUADRATIC
#define SGL_DEBUG_QUADRATIC_STOPPING
#endif

#ifdef SGL_DEBUG_BLOCK_ACTIVE
#define SGL_SHOW_HEAVY_DEBUG_WARNING
#endif
#ifdef SGL_DEBUG_GB
#define SGL_SHOW_HEAVY_DEBUG_WARNING
#endif

#endif /* CONFIG_H_ */
