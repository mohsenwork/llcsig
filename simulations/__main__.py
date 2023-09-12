import logging
import sys
from multiprocessing import Pool
from .multiprocess_simulations import Multiprocessing
from .optuna_custome import OptunaCustome
from .objective_1 import Objective1
from .objective_2 import Objective2
from .objective_3 import Objective3
from .objective_4 import Objective4
from .objective_5 import Objective5
from .objective_6 import Objective6
from .objective_7 import Objective7

if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                        handlers=[
                            logging.FileHandler('logfile.log'),
                            logging.StreamHandler(sys.stdout)
                        ],
                        datefmt='%H:%M:%S',
                        level=logging.INFO)
    logger = logging.getLogger(__name__)

    flag = int(sys.argv[1])

    if flag == 0:
        Multiprocessing.run(*sys.argv[2:])

    # LLC-NF using ROC AUC
    if flag == 1:
        OptunaCustome.run(Objective1, *sys.argv[2:])

    # LLC-NF using accuracy
    if flag == 2:
        OptunaCustome.run(Objective2, *sys.argv[2:])

    # LLC-F using ROC AUC
    if flag == 3:
        OptunaCustome.run(Objective3, *sys.argv[2:])
    
    # LLC-F using accuracy
    if flag == 4:
        OptunaCustome.run(Objective4, *sys.argv[2:])
    
    # ASP using ROC AUC
    if flag == 5:
        OptunaCustome.run(Objective5, *sys.argv[2:])

    # ASP using accuracy
    if flag == 6:
        OptunaCustome.run(Objective6, *sys.argv[2:])
    
    # ASP narrow using accuracy
    if flag == 7:
        OptunaCustome.run(Objective7, *sys.argv[2:])

    logger.info('Done')