#include <string.h>
#include "pdc_cont.h"
#include "pdc_malloc.h"
#include "pdc_prop_pkg.h"
#include "pdc_interface.h"

static perr_t PDCcont__close(struct PDC_cont_info *cp);

/* PDC container ID class */
static const PDCID_class_t PDC_CONT_CLS[1] = {{
    PDC_CONT,                           /* ID class value */
    0,                                  /* Class flags */
    0,                                  /* # of reserved IDs for class */
    (PDC_free_t)PDCcont__close          /* Callback routine for closing objects of this class */
}};

perr_t PDCcont_init(PDC_CLASS_t *pc)
{
    perr_t ret_value = SUCCEED;         /* Return value */

    FUNC_ENTER(NULL);

    /* Initialize the atom group for the container IDs */
    if(PDC_register_type(PDC_CONT_CLS, pc) < 0)
        PGOTO_ERROR(FAIL, "unable to initialize container interface");

done:
    FUNC_LEAVE(ret_value); 
} /* end PDCcont_init() */

pdcid_t PDCcont_create(pdcid_t pdc, const char *cont_name, pdcid_t cont_create_prop)
{
    pdcid_t ret_value = SUCCEED;
    struct PDC_cont_info *p = NULL;
    pdcid_t new_id;

    FUNC_ENTER(NULL);

    p = PDC_MALLOC(struct PDC_cont_info);
    if(!p)
        PGOTO_ERROR(FAIL,"PDC container memory allocation failed\n");
    p->name = strdup(cont_name);
    p->pdc = pdc;
    p->cont_prop = cont_create_prop;
//    p->cont_life = PDC_PERSIST;
    new_id = PDC_id_register(PDC_CONT, p, pdc);
    ret_value = new_id;
    
done:
    FUNC_LEAVE(ret_value);
} /* end of PDCcont_create() */

perr_t PDC_cont_list_null(pdcid_t pdc)
{
    perr_t ret_value = SUCCEED;   /* Return value */
    int nelemts;
    
    FUNC_ENTER(NULL);
    // list is not empty
    nelemts = PDC_id_list_null(PDC_CONT, pdc);
    if(PDC_id_list_null(PDC_CONT, pdc) > 0) {
        /* printf("%d element(s) in the container list will be automatically closed by PDC_close()\n", nelemts); */
        if(PDC_id_list_clear(PDC_CONT, pdc) < 0)
            PGOTO_ERROR(FAIL, "fail to clear container list");
    }

done:
    FUNC_LEAVE(ret_value);
}

static perr_t PDCcont__close(struct PDC_cont_info *cp)
{
    perr_t ret_value = SUCCEED;         /* Return value */

    FUNC_ENTER(NULL);

    cp = PDC_FREE(struct PDC_cont_info, cp);
    
    FUNC_LEAVE(ret_value);
} /* end of PDCcont__close() */

perr_t PDCcont_close(pdcid_t id, pdcid_t pdc)
{
    perr_t ret_value = SUCCEED;   /* Return value */

    FUNC_ENTER(NULL);

    /* When the reference count reaches zero the resources are freed */
    if(PDC_dec_ref(id, pdc) < 0)
        PGOTO_ERROR(FAIL, "container: problem of freeing id");
    
done:
    FUNC_LEAVE(ret_value);
} /* end of PDCcont_close() */

perr_t PDCcont_end(pdcid_t pdc)
{
    perr_t ret_value = SUCCEED;         /* Return value */

    FUNC_ENTER(NULL);

    if(PDC_destroy_type(PDC_CONT, pdc) < 0)
        PGOTO_ERROR(FAIL, "unable to destroy container interface");
    
done:
    FUNC_LEAVE(ret_value);
} /* end of PDCcont_end() */

pdcid_t PDCcont_open(pdcid_t pdc_id, const char *cont_name)
{
    pdcid_t ret_value = SUCCEED;
    pdcid_t cont_id;

    FUNC_ENTER(NULL);

    // should wait for response from server 
    // look up in the list for now
    cont_id = PDC_find_byname(PDC_CONT, cont_name, pdc_id);
    if(cont_id <= 0)
        PGOTO_ERROR(FAIL, "cannot locate container");
    pdc_inc_ref(cont_id, pdc_id);
    ret_value = cont_id;
    
done:
    FUNC_LEAVE(ret_value);
} /* end of PDCcont_open() */

cont_handle *PDCcont_iter_start(pdcid_t pdc_id)
{
    cont_handle *ret_value = NULL;
    cont_handle *conthl = NULL;
    PDC_CLASS_t *pc;
    struct PDC_id_type *type_ptr;

    FUNC_ENTER(NULL);

    pc = (PDC_CLASS_t *)pdc_id;
    type_ptr  = (pc->PDC_id_type_list_g)[PDC_CONT];
    if(type_ptr == NULL) 
        PGOTO_ERROR(NULL, "container list is empty");
    conthl = (&type_ptr->ids)->head;
    ret_value = conthl;
    
done:
    FUNC_LEAVE(ret_value);
} /* end of PDCcont_iter_start() */

pbool_t PDCcont_iter_null(cont_handle *chandle)
{
    pbool_t ret_value = FALSE;
    
    FUNC_ENTER(NULL);
    
    if(chandle == NULL)
        ret_value = TRUE;
    
    FUNC_LEAVE(ret_value); 
} /* end of PDCcont_iter_null() */

cont_handle *PDCcont_iter_next(cont_handle *chandle)
{
    cont_handle *ret_value = NULL;
    cont_handle *next = NULL;

    FUNC_ENTER(NULL);

    if(chandle == NULL)
        PGOTO_ERROR(NULL, "no next container");
    next = PDC_LIST_NEXT(chandle, entry); 
    ret_value = next;
    
done:
    FUNC_LEAVE(ret_value);
} /* end of PDCcont_iter_next() */

struct PDC_cont_info *PDCcont_iter_get_info(cont_handle *chandle)
{
    struct PDC_cont_info *ret_value = NULL;
    struct PDC_cont_info *info = NULL;

    FUNC_ENTER(NULL);

    info = (struct PDC_cont_info *)(chandle->obj_ptr);
    if(info == NULL)
        PGOTO_ERROR(NULL, "PDC container info memory allocation failed");
    
    ret_value = info;
    
done:
    FUNC_LEAVE(ret_value);
} /* end of PDCcont_iter_get_info() */

perr_t PDCcont_persist(pdcid_t cont_id, pdcid_t pdc)
{
    perr_t ret_value = SUCCEED;         /* Return value */
    PDC_CLASS_t *pc;
    pdcid_t propid;
    struct PDC_id_info *info;
    struct PDC_id_info *prop;
    
    FUNC_ENTER(NULL);
    
    pc = (PDC_CLASS_t *)pdc;
    info = PDC_find_id(cont_id, pc);
    if(info == NULL)
        PGOTO_ERROR(FAIL, "cannot locate container ID");
    propid = ((struct PDC_cont_info *)(info->obj_ptr))->cont_prop;
    prop = PDC_find_id(propid, pc);
    if(prop == NULL)
        PGOTO_ERROR(FAIL, "cannot container property");
    ((struct PDC_cont_prop *)(prop->obj_ptr))->cont_life = PDC_PERSIST;
    
done:
    FUNC_LEAVE(ret_value);
} /* end of PDCcont_persist() */

perr_t PDCprop_set_cont_lifetime(pdcid_t cont_prop, PDC_lifetime cont_lifetime, pdcid_t pdc)
{
    perr_t ret_value = SUCCEED;         /* Return value */
    PDC_CLASS_t *pc;
    struct PDC_id_info *info;
    
    FUNC_ENTER(NULL);
    
    pc = (PDC_CLASS_t *)pdc;
    info = PDC_find_id(cont_prop, pc);
    if(info == NULL)
        PGOTO_ERROR(FAIL, "cannot locate container property ID");
    ((struct PDC_cont_prop *)(info->obj_ptr))->cont_life = cont_lifetime;
    
done:
    FUNC_LEAVE(ret_value);
} /* end of PDCprop_set_cont_lifetime() */
